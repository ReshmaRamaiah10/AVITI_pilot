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

#ifndef REFSEQ_H_
#define REFSEQ_H_

#include <cctype>
#include <cstdint>
#include <fstream>
#include <string>

#include "utils.h"
#include "my_assert.h"

class RefSeq {
public:
	RefSeq();
	RefSeq(const std::string& name, const std::string& seq);

	void set(const std::string& name, const std::string& seq) {
		this->name = name; this->seq = seq;
		len = seq.length(); assert(len > 0);
	}
	
	bool read(std::ifstream& fin);
	void write(std::ofstream& fout);

	int getLen() const { return len; }
	const std::string& getName() const { return name; }
	const std::string& getSeq() const { return seq; }

	char baseAt(char dir, int pos) const {
		general_assert(pos >= 0 && pos < len, "At RefSeq.baseAt. Requested position is out of boundary: pos = " + itos(pos) + ", transcript = " + name + ", length = " + itos(len) + ".");
		return (dir == '+' ? seq[pos] : base2rbase[seq[len - pos - 1]]);
	}

	int baseCodeAt(char dir, int pos) const {
		general_assert(pos >= 0 && pos < len, "At RefSeq.baseCodeAt. Requested position is out of boundary: pos = " + itos(pos) + ", transcript = " + name + ", length = " + itos(len) + ".");
		return (dir == '+' ? base2code[seq[pos]] : rbase2code[seq[len - pos - 1]]);
	}
	
private:
	int len; // the transcript length
	std::string name; // the tag
	std::string seq; // the sequence, in the forward strand	
};

#endif
