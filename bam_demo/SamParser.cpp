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

#include <cstdio>
#include <cstring>
#include <string>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "my_assert.h"
#include "SamHeaderText.hpp"
#include "SamParser.hpp"

SamParser::SamParser(const char* inpF, htsThreadPool* p) {
	sam_in = sam_open(inpF, "r");
	general_assert(sam_in != 0, "Cannot open " + cstrtos(inpF) + "! It may not exist.");

	header = sam_hdr_read(sam_in);
	general_assert(header != 0, "Fail to parse the header!");
	ht = new SamHeaderText(header);

	if (p != NULL) hts_set_opt(sam_in, HTS_OPT_THREAD_POOL, p);
}

SamParser::~SamParser() {
	delete ht;
	bam_hdr_destroy(header);
	sam_close(sam_in);
}
