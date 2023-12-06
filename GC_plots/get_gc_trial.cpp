#include <iostream>
#include <fstream>
#include <vector>
#include <zlib.h>
#include <cstring>

using namespace std;

int main() {
    // Initialize counters
    vector<int> gc_counts(100, 0); // Counter array for percent of GCs in 0-1%, 1-2%, ..., 99-100%.

    // Open R2 read files for lane 1 and lane 2
    gzFile r2_lane1 = gzopen("/home/ramaiahr/AVITI_fastqs/Ch10_scRNA_AVITI/Ch10_scRNA_AVITI_S6_L001_R2_001.fastq.gz", "r");
    gzFile r2_lane2 = gzopen("/home/ramaiahr/AVITI_fastqs/Ch10_scRNA_AVITI/Ch10_scRNA_AVITI_S6_L002_R2_001.fastq.gz", "r");

    // Loop over each read in lane 1
    string line;
    line.resize(500); // Set the size of the line string to the maximum length of a line in the FASTQ file plus one
    int line_num = 0; // Keeps track of which line of the FASTQ file we're on
    while (gzgets(r2_lane1, &line[0], line.size())) {
        line_num++;
        if (line_num % 4 == 2) { // Only process the sequence line
            int gc_count = 0;
            int line_length = strlen(&line[0]); // Get the actual length of the line
            for (int i = 0; i < line_length; i++) {
                if (line[i] == 'G' || line[i] == 'C') {
                    gc_count++;
                }
            }
            int percent_gc = gc_count * 100 / line_length;
            gc_counts[percent_gc]++;
        }
    }
    // Loop over each read in lane 2
    line_num = 0;
    while (gzgets(r2_lane2, &line[0], line.size())) {
        line_num++;
        if (line_num % 4 == 2) { // Only process the sequence line
            int gc_count = 0;
            int line_length = strlen(&line[0]); // Get the actual length of the line
            for (int i = 0; i < line_length; i++) {
                if (line[i] == 'G' || line[i] == 'C') {
                    gc_count++;
                }
            }
            int percent_gc = gc_count * 100 / line_length;
            gc_counts[percent_gc]++;
        }
    }

    // Output results to text file
    ofstream output("Ch10_scRNA_AVITI_gc_counts.txt");
    for (int i = 0; i < 100; i++) {
        output << i << "-" << (i+1) << "%\t" << gc_counts[i] << endl;
    }
    output.close();

    // Generate curves based on counts (not included here)

    // Close files
    gzclose(r2_lane1);
    gzclose(r2_lane2);

    return 0;
}
