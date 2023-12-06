#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>

#include "gzip_utils.hpp"

// Command example: ./get_qc 200 R1_L1.fastq.gz R1_L2.fastq.gz ./result/qc.csv

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cout << "Missing arguments!" << std::endl;
    return 1;
  }

  int read_length = std::stoi(std::string(argv[1]));
  if (read_length > 0)
    std::cout << "Only consider the first " << read_length << " positions." << std::endl;
  else
    read_length = std::numeric_limits<int>::max();


  Read read;
  // Quality scores: 0 - 40 for Illumina and 0 - 50 for AVITI
  int max_quality_score = 50;
  std::vector<std::vector<int> > qual_count;
  int read_count;

  for (int i = 2; i < argc - 1; ++i) {
    std::cout << i - 1 << ". Processing " << argv[i] << std::endl;
    iGZipFile gzip_in(argv[i]);

    // Iterate over lines.
    int len;
    read_count = 0;

    while (gzip_in.next(read)) {
      len = std::min(read_length, int(read.seq.length()));
      for (int j = 0; j < len; ++j) {
        if (qual_count.size() < j + 1)
          qual_count.push_back(std::vector<int>(max_quality_score + 1, 0));
        ++qual_count[j][read.qual[j] - 0x21];
      }
      ++read_count;
    }

    std::cout << "Total reads = " << read_count << "." << std::endl;
  }

  // Write statistics into file.
  std::ofstream fout;
  fout.open(argv[argc - 1]);
  for (int i = 0; i < qual_count.size(); ++i) {
    for (int j = 0; j <= max_quality_score; ++j) {
      fout << qual_count[i][j];
      if (j < max_quality_score)
        fout << ",";
      else
        fout << std::endl;
    }
  }
  fout.close();


  return 0;
}
