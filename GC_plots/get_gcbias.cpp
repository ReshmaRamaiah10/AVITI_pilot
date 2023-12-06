#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "gzip_utils.hpp"

// Command example: ./get_gc_stats sample_lane1.fastq.gz sample_lane2.fastq.gz ./result/gc_bias.csv

  int main(int argc, char** argv) {
  if (argc < 4) {
    std::cout << "Missing arguments!" << std::endl;
    return 1;
  }

  std::vector<int> gc_counts(100, 0); // counter array for percent of GCs in 0-1%, 1-2%, ..., 99-100%

  // Iterate lines
  std::string line;
  int read_length;
  Read read;
  for (int i = 1; i < argc - 1; ++i) {
    iGZipFile gzip_in(argv[i]);
    std::cout << " processing " << argv[i] << std::endl;
    while (gzip_in.next(read)) {
      read_length = read.seq.length();
      int gc_count = 0;
      for (char c : read.seq) {
        if (c == 'G' || c == 'C') {
          ++gc_count;
        }
      }
      int gc_percent = gc_count * 100 / read_length;
      ++gc_counts[gc_percent];
    }
  }

  // Write statistics into file.
  std::ofstream fout;
  fout.open(argv[argc - 1]);
  for (int i = 0; i < gc_counts.size(); ++i) {
    fout << i << "-" << i+1 << "%," << gc_counts[i] << std::endl;
  }
  fout.close();

  return 0;
}
