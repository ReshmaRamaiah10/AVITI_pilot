/* Copyright (c) 2017
   Bo Li (The Broad Institute of MIT and Harvard)
   libo@broadinstitute.org

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.   

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*/

#ifndef BAMALIGNMENT_H_
#define BAMALIGNMENT_H_

#include <cmath>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <string>
#include <algorithm>

#include "htslib/sam.h"

#include "my_assert.h"
#include "CIGARstring.hpp"
#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "MDstring.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"

class BamAlignment {
public:
	BamAlignment() : is_paired(false), is_primary(false), is_aligned(-1),  b(NULL), b2(NULL) {
	}

	BamAlignment(const BamAlignment& o) : b(NULL), b2(NULL) {
		is_paired = o.is_paired;
		is_primary = o.is_primary;
		is_aligned = o.is_aligned;
		if (o.b != NULL) b = bam_dup1(o.b);
		if (o.b2 != NULL) b2 = bam_dup1(o.b2);
	}

	~BamAlignment() {
		if (b != NULL) bam_destroy1(b);
		if (b2 != NULL) bam_destroy1(b2);
	}
	
	/*
		@param   in   input SamParser 
	 */
	bool read(SamParser* in);

	/*
		@param   out     output BamWriter
	 */
	bool write(BamWriter* out);
	
	// overall stats
	
	bool isPaired() const { return is_paired; }
	
	// -1: nothing is loaded, 0: not aligned, 1: first mate is aligned (including SE reads), 2: second mate is aligned, 3: fully aligned 
	int isAligned() const { return is_aligned; }
	
	bool isPrimary() const { return is_primary; }

	int getInsertSize() const {
		assert(is_aligned == 3);
		return abs(b->core.isize); 
	}

	// 0-based, directly load from BAM/SAM file
	int get_tid() const { 
		assert(is_aligned > 0); 
		return (is_aligned & 1) ? b->core.tid : b2->core.tid;
	}

	// 1-based, used for indexing Transcripts and Refs 
	int getTid() const {
		return get_tid() + 1;
	}
	
	char getDir() const { 
		assert(is_aligned > 0);
		if (is_aligned & 1) return !bam_is_rev(b) ? '+' : '-';
		return bam_is_rev(b2) ? '+' : '-';
	}

	char getMateDir(int mate = 1) const {
		assert((mate == 1 && bam_is_mapped(b)) || (mate == 2 && is_paired && bam_is_mapped(b2)));
		if (mate == 1) return bam_is_rev(b) ? '-' : '+';
		else return bam_is_rev(b2) ? '-' : '+';
	}

	/*
	  @param     fragment_length     The average fragment length, 0 means no fragment length is provided
	  @return    If fragment_length == 0, return the leftmost position of two mates. Otherwise, return the leftmost position calculated with fragment length.
	             If the calculated position < 0, set it to 0
	  @comment   bam_endpos gives the right-most position + 1 
	 */
	int getLeftMostPos(int fragment_length = 0) const {
		assert(is_aligned > 0);
		if (is_aligned == 3) return std::min(b->core.pos, b2->core.pos);
		if (fragment_length <= 0) return (is_aligned & 1) ? b->core.pos : b2->core.pos;

		if (is_aligned & 1) return bam_is_rev(b) ? std::max(0, int(bam_endpos(b)) - fragment_length) : b->core.pos;
		else return bam_is_rev(b2) ? std::max(0, int(bam_endpos(b2)) - fragment_length) : b2->core.pos;
	}
		
	int getMapQ() const { return (is_aligned & 1) ? b->core.qual : b2->core.qual; }

	void setMapQ(int MapQ) { 
		b->core.qual = (bam_is_mapped(b) ? MapQ : 255);
		if (is_paired) b2->core.qual = (bam_is_mapped(b2) ? MapQ : 255);
	}

	// for mates
	int getPos(int mate) const {
		assert(is_aligned > 0 && (mate == 1 || (is_paired && mate == 2)));
		return (mate == 1 ? b->core.pos : b2->core.pos);
	}
	
	/*
	  @param   mate   which mate (1 or 2)
	  @param   target_len   the length of reference sequence
	  @comment: This function returns the smallest position of the mate from its strand
	 */
	int getDirPos(int mate, int target_len) const {
		assert(is_aligned > 0 && (mate == 1 || (is_paired && mate == 2)));
		if (mate == 1) return bam_is_rev(b) ? target_len - bam_endpos(b) : b->core.pos;
		else return bam_is_rev(b2) ? target_len - bam_endpos(b2) : b2->core.pos;
	}

	// mate = 0, the name of the read; 1, mate 1; 2, mate 2
	const char* getName(int mate = 0) const { 
		return bam_get_qname(mate < 2 ? b : b2); 
	}

	// length of the query sequence
	int getSeqLength(int mate = 1) const { 
		assert(mate == 1 || (is_paired && mate == 2));
		assert(((mate == 1) && (b->core.l_qseq > 0)) || ((mate == 2) && (b2->core.l_qseq > 0)));
		return mate == 1 ? b->core.l_qseq : b2->core.l_qseq;
	} 
	
	// length of the aligned reference segment	
	int getAlignedLength(int mate = 1) const {
		assert(mate == 1 || (is_paired && mate == 2));
		return mate == 1 ? bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)) : bam_cigar2rlen(b2->core.n_cigar, bam_get_cigar(b2));
	}
	
	bool getCIGAR(CIGARstring& ci, int mate = 1) {
		assert(mate == 1 || (is_paired && mate == 2));
		if (mate == 1) ci.setUp(bam_get_cigar(b), b->core.n_cigar, is_ori(b));
		else ci.setUp(bam_get_cigar(b2), b2->core.n_cigar, is_ori(b2));
		return true;
	}
	
	bool getSEQ(SEQstring& si, int mate = 1) {
		assert(mate == 1 || (is_paired && mate == 2));
		assert(((mate == 1) && (b->core.l_qseq > 0)) || ((mate == 2) && (b2->core.l_qseq > 0)));
		if (mate == 1) si.setUp(bam_get_seq(b), b->core.l_qseq, is_ori(b));
		else si.setUp(bam_get_seq(b2), b2->core.l_qseq, is_ori(b2));
		return true;
	}
	
	bool getQUAL(QUALstring& qi, int mate = 1) {
		assert(mate == 1 || (is_paired && mate == 2));
		assert(((mate == 1) && (b->core.l_qseq > 0)) || ((mate == 2) && (b2->core.l_qseq > 0)));
		if (mate == 1) qi.setUp(bam_get_qual(b), b->core.l_qseq, is_ori(b));
		else qi.setUp(bam_get_qual(b2), b2->core.l_qseq, is_ori(b2));
		return true;
	}
	
	bool getMD(MDstring& mdstr, int mate = 1) {
		uint8_t* p;    
		assert(mate == 1 || (is_paired && mate == 2));
		p = bam_aux_get((mate == 1 ? b : b2), "MD");
		assert(p != NULL && *p++ == 'Z');
		mdstr.setUp((char*)p);
		return true;
	}

	// optional fields
	char tag2A(const uint8_t* p) const { return bam_aux2A(p); }
	int tag2i(const uint8_t* p) const { return bam_aux2i(p); }
	unsigned int tag2u(const uint8_t* p) const { return bam_aux2i(p); }	
	double tag2f(const uint8_t* p) const { return bam_aux2f(p); }
	char* tag2Z(const uint8_t* p) const { return bam_aux2Z(p); }

	// type can be either 'f': float, 'i': signed int, 'u': unsigned int, 'Z': char array, 'A': char, or 'B' : byte array
	bool findTag(const char tag[2], uint8_t*& p, char& type, int mate = 1) {
		assert(mate == 1 || (is_paired && mate == 2));
		p = bam_aux_get((mate == 1 ? b : b2), tag);
		if (p == NULL) return false;
		type = *p;
		if (type == 'c' || type == 's' || type == 'i') type = 'i';
		else if (type == 'C' || type == 'S' || type == 'I') type = 'u';
		else if (type == 'd' || type == 'f') type = 'f';
		else assert(type == 'A' || type == 'Z' || type == 'H' || type == 'B');
		return true;
	}

	void insertTag(const char tag[2], char type, int len, const uint8_t* data, int mate = 0) {
		uint8_t *p_tag = NULL;
		int old_len, delta, rest_len;

		if (mate != 2) {
			p_tag = bam_aux_get(b, tag);
			if (p_tag == NULL) bam_aux_append(b, tag, type, len, data);
			else {
				assert(*p_tag == type);
				old_len = bam_aux_type2size(*p_tag, p_tag);
				delta = p_tag - b->data;
				rest_len = b->l_data - (delta + 1 + old_len);
				b->l_data += len - old_len; expand_data_size(b);
				p_tag = b->data + delta; // reassign since expand_data_size might realloc b->data
				if (old_len != len && rest_len > 0) memmove(p_tag + 1 + len, p_tag + 1 + old_len, rest_len);
				memcpy(p_tag + 1, data, len);
			}			
		}

		if (is_paired && mate != 1) {
			p_tag = bam_aux_get(b2, tag);
			if (p_tag == NULL) bam_aux_append(b2, tag, type, len, data);
			else {
				assert(*p_tag == type);
				old_len = bam_aux_type2size(*p_tag, p_tag);
				delta = p_tag - b2->data;
				rest_len = b2->l_data - (delta + 1 + old_len);
				b2->l_data += len - old_len; expand_data_size(b2);
				p_tag = b2->data + delta; // reassign since expand_data_size might realloc b->data
				if (old_len != len && rest_len > 0) memmove(p_tag + 1 + len, p_tag + 1 + old_len, rest_len);
				memcpy(p_tag + 1, data, len);
			}
		}
	}

	void removeTag(const char tag[2], int mate = 0) {
		uint8_t *p_tag = NULL;

		if (mate != 2 && (p_tag = bam_aux_get(b, tag)) != NULL) bam_aux_del(b, p_tag);
		if (is_paired && mate != 1 && (p_tag = bam_aux_get(b2, tag)) != NULL) bam_aux_del(b2, p_tag);
	}

private:
	static const uint8_t rnt_table[16];

	bool is_paired, is_primary;
	char is_aligned; // 2 bits, from right to left, the first bit represents first mate and the second bit represents the second mate
                   	 // Thus, 0, unalignable; 1, only first mate; 2, only second mate; 3, both mates
	
	bam1_t *b, *b2;
	
	bool bam_is_paired(const bam1_t* b) const { return (b->core.flag & BAM_FPAIRED); }
	bool bam_is_proper(const bam1_t* b) const { return (b->core.flag & BAM_FPROPER_PAIR); }
	bool bam_is_mapped(const bam1_t* b) const { return !(b->core.flag & BAM_FUNMAP); }
	bool bam_is_unmapped(const bam1_t* b) const { return (b->core.flag & BAM_FUNMAP); }
	bool bam_is_read1(const bam1_t* b) const { return (b->core.flag & BAM_FREAD1); }
	bool bam_is_read2(const bam1_t* b) const { return (b->core.flag & BAM_FREAD2); }
	bool bam_is_primary(const bam1_t* b) const { return !(b->core.flag & BAM_FSECONDARY); }
	
	int bam_aux_type2size(char x, const uint8_t* p = NULL) {
		if (x == 'C' || x == 'c' || x == 'A') return 1;
		else if (x == 'S' || x == 's') return 2;
		else if (x == 'I' || x == 'i' || x == 'f') return 4;
		else if (x == 'd') return 8;
		else if (x == 'Z' || x == 'H') { assert(p != NULL); return strlen((char*)(p + 1)) + 1; }
		else {
			assert(x == 'B' && p != NULL);
			return 1 + 4 + bam_auxB_len(p) * bam_aux_type2size(*(p + 1));
		}
	}

	// if the sequence and qual are original
	bool is_ori(bam1_t* b) {
		return ((bam_is_unmapped(b) || !bam_is_rev(b)) ? true : false);
	}

	// Caution: this function may change b->data's adddress!  
	void expand_data_size(bam1_t* b) {
		if (b->m_data < b->l_data) {
			b->m_data = b->l_data;
			kroundup32(b->m_data);
			b->data = (uint8_t*)realloc(b->data, b->m_data);
		}
	}
};

#endif
