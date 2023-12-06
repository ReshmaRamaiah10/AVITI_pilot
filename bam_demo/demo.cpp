#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cassert>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "htslib/hts.h"

#include "Refs.hpp"
#include "SamParser.hpp"
#include "BamAlignment.hpp"
using namespace std;

const int max_qlen=50;

bool verbose = true;

int num_threads;
htsThreadPool p = {NULL, 0};

Refs refs;
SamParser *parser;
BamAlignment ba;

std::vector<long> mis, all;
long insertion_length, deletion_length, processed_length;

// Make sure the sequence name in the genome FASTA is identical to the BAM header
void check_seqname() {
	const bam_hdr_t* header = parser->getHeader();
	assert(header->n_targets == refs.getM());
	for (int i = 0; i < header->n_targets; ++i)
		assert(!strcmp(header->target_name[i], refs.getRef(i).getName().c_str()));
}


inline void collect_data(char dir, int pos, const RefSeq& refseq, const CIGARstring& cigar, const SEQstring& seq, const QUALstring& qual) {
	int len = cigar.getLen();
	char opchr;
	int oplen;
	int refpos = pos, readpos = 0;
	int ref_base, read_base, qual_score;

	for (int i = 0; i < len; ++i) {
		opchr = cigar.opchrAt(i);
		oplen = cigar.oplenAt(i);

		if (opchr == 'M' || opchr == '=' || opchr == 'X') {
			for (int j = 0; j < oplen; ++j) {
				ref_base = refseq.baseCodeAt(dir, refpos);
				read_base = seq.baseCodeAt(readpos);
				qual_score = qual.qualAt(readpos);
				++all[qual_score];
				mis[qual_score] += (ref_base != read_base);
				++refpos; ++readpos;
			}
		}
		else if (opchr == 'I' || opchr == 'S') {
			insertion_length += oplen;
			readpos += oplen;
		}
		else if (opchr == 'D' || opchr == 'N') {
			deletion_length += oplen;
			refpos += oplen;
		}
	}
}


int main(int argc, char* argv[]) {
	if (argc != 5) {
		printf("Usage: demo input.fa input.bam num_threads output_txt\n");
		exit(-1);
	}

	num_threads = atoi(argv[3]);
	refs.readFrom(argv[1]);
	if (num_threads > 1) assert(p.pool = hts_tpool_init(num_threads));
	parser = new SamParser(argv[2], num_threads > 1 ? &p : NULL);

	check_seqname();

	mis.resize(max_qlen, 0.0);
	all.resize(max_qlen, 0.0);
	insertion_length = 0;
	deletion_length = 0;
	processed_length = 0;

	long cnt = 0;
	char dir;
	int pos, tid;
	CIGARstring ci;
	SEQstring si;
	QUALstring qi;

	while (ba.read(parser)) {
		if ((ba.isAligned() & 1) && ba.isPrimary() && (ba.getMapQ() == 255)) {
			dir = ba.getDir();
			tid = ba.get_tid();
			pos = ba.getDirPos(1, refs.getRef(tid).getLen());
			ba.getCIGAR(ci);
			ba.getSEQ(si);
			ba.getQUAL(qi);
			collect_data(dir, pos, refs.getRef(tid), ci, si, qi);
			processed_length += ba.getSeqLength();
		}
		++cnt;
		if (cnt % 1000000==0) {
			printf("Processed %ld lines\n", cnt);
			//break;
		}
	}

	for (int i = 0; i < max_qlen; ++i) {
		printf("%d => %ld\n", i, all[i]);
	}

	printf("Insertion %ld, Deletion %ld, Total processed %ld\n", insertion_length, deletion_length, processed_length);

	FILE *fo = fopen(argv[4], "w");
	for (int i = 0; i < max_qlen; ++i) {
		if (all[i] > 0) {
			double p = mis[i] * 1.0 / all[i];
			fprintf(fo, "%d,%d,%ld,%ld\n", i, int(-10 * log10(p) + 0.5), mis[i], all[i]);
		}
	}
	fclose(fo);

	delete parser;
	if (num_threads > 1) hts_tpool_destroy(p.pool);

	return 0;
}
