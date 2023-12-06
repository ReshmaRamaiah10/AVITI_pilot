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
#include <cassert>
#include <string>
#include <fstream>
#include <vector>

#include "my_assert.h"
#include "RefSeq.hpp"
#include "Refs.hpp"

Refs::Refs() {
	M = 0;
	seqs.clear();
}

void Refs::readFrom(char* inpF) {
	std::ifstream fin(inpF);
	std::string line, chrname, seqname, gseq;
	int seqlen, len;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(inpF) + "! It may not exist.");

	M = 0;
	seqs.clear();	
	getline(fin, line);
	while ((fin) && (line[0] == '>')) {
		std::istringstream strin(line.substr(1));
		strin>> chrname>> seqname;
		
		gseq = ""; seqlen = 0;
		while((getline(fin, line)) && (line[0] != '>')) {
		  len = line.length();
		  seqlen += len;
		  gseq += line;
		}
		assert(seqlen > 0);

		seqs.emplace_back(chrname, gseq);
		++M;
	}

	fin.close();
	assert(M == (int)seqs.size());

	if (verbose) printf("%s is loaded.\n", inpF);
}

void Refs::writeTo(char* outF) {
	std::ofstream fout(outF);
	for (int i = 1; i <= M; ++i) 
		seqs[i].write(fout);
	fout.close();
	if (verbose) printf("%s is generated.\n", outF);
}
