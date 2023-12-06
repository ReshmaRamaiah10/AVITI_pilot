# run fastqc on a fastq.gz files
echo "Running fastqc for Ch7_scRNA_AVITI_S5_L001_R1_001.fastq.gz"
fastqc AVITI_fastq_files/Ch7_scRNA_AVITI/Ch7_scRNA_AVITI_S5_L001_R1_001.fastq.gz -o AVITI_fastqc_files/Ch7_scRNA_AVITI
# OR
# run fastqc on all fastq.gz files in a directory
echo "Running fastqc for Ch7_scRNA_AVITI"
fastqc AVITI_fastq_files/Ch7_scRNA_AVITI/* -o AVITI_fastqc_files/Ch7_scRNA_AVITI

# Once all fastqc files are generated for a sample run multiqc to combine all fastqc results
multiqc AVITI_fastqc_files/Ch7_scRNA_AVITI/* -o multiqc_results/Ch7_scRNA_AVITI
