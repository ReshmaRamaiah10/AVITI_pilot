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

#include <cmath>
#include <cstdio>
#include <cctype>
#include <cstring>
#include <cstdint>
#include <cassert>
#include <string>
#include <algorithm>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "utils.h"
#include "my_assert.h"
#include "SamParser.hpp"
#include "BamWriter.hpp"
#include "BamAlignment.hpp"

const uint8_t BamAlignment::rnt_table[16] = {0, 8, 4, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 15};

bool BamAlignment::read(SamParser *in) {
	is_aligned = -1;

	if (b == NULL) b = bam_init1();
	if (!in->read(b)) return false;

	is_paired = bam_is_paired(b);
	assert(!is_paired);
	// read the second mate
	if (is_paired) { 
		if (b2 == NULL) b2 = bam_init1();

		general_assert(in->read(b2) && bam_is_paired(b2), "Fail to read the other mate for a paired-end alignment!");
		general_assert(((b->core.flag & 0x00C0) == 0x0040 && (b2->core.flag & 0x00C0) == 0x0080) || 
			 ((b->core.flag & 0x00C0) == 0x0080 && (b2->core.flag & 0x00C0) == 0x0040), 
			 "Cannot detect both mates of a paired-end alignment!");
		
		if (bam_is_read2(b)) { bam1_t* tmp = b; b = b2; b2 = tmp; }
	}

	// calculate is_aligned
	is_aligned = bam_is_mapped(b);
	if (is_paired) is_aligned |= ((char)bam_is_mapped(b2) << 1);	
	is_primary = (is_aligned == 2) ? bam_is_primary(b2): bam_is_primary(b);

	// The following four statements are grouped together
	int tmp_len = 0, tmp_len2 = 0;
	tmp_len = b->core.l_qseq;
	if (is_paired) tmp_len2 = b2->core.l_qseq;
	assert(!(is_aligned & 1) || tmp_len == bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)));
	assert(!(is_aligned & 2) || tmp_len2 == bam_cigar2qlen(b2->core.n_cigar, bam_get_cigar(b2)));

	return true;
}

bool BamAlignment::write(BamWriter *out) {
	assert(is_aligned >= 0 && b != NULL && (!is_paired || b2 != NULL));
	out->write(b);
	if (is_paired) out->write(b2);

	return true;
}
