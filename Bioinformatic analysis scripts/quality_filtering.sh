#! /usr/bin/env bash
#
# Author: Sarah François
# Date of last update: 2019-10-07
#
# Usage: ./quality_filtering.sh SAMPLE_LIST
#
# EXAMPLE: ./quality_filtering.sh sample_list.txt
#

# SET VARIABLES
SAMPLELIST=$1

# CREATE DIRECTORY FOR RESULTS
mkdir filtered_reads

# DATA CLEANING
cat $SAMPLELIST | while read line
do
# Add sample name to reads IDs
  unpigz -p 24 $line.1.fastq.gz # Uncompress R1 file
  unpigz -p 24 $line.2.fastq.gz # Uncompress R2 file
  paste - - - - < $line.1.fastq > filtered_reads/$line.1.tab # Convert fastq to tab
  paste - - - - < $line.2.fastq > filtered_reads/$line.2.tab # Convert fastq to tab
  awk -F "\t" -v mid=$line '{print $0"\t"mid}' filtered_reads/$line.1.tab >> filtered_reads/$line.1.name.tab # Add sample name to reads ids
  awk -F "\t" -v mid=$line '{print $0"\t"mid}' filtered_reads/$line.2.tab >> filtered_reads/$line.2.name.tab # Add sample name to reads ids
    rm filtered_reads/$line.1.tab
    rm filtered_reads/$line.2.tab
  awk -F "\t" '{print $1"\t"$5"\n"$2"\n"$3"\n"$4}' filtered_reads/$line.1.name.tab > filtered_reads/$line.1.name.fastq # Convert tab to fastq
  awk -F "\t" '{print $1"\t"$5"\n"$2"\n"$3"\n"$4}' filtered_reads/$line.2.name.tab > filtered_reads/$line.2.name.fastq # Convert tab to fastq
    rm filtered_reads/$line.1.name.tab
    rm filtered_reads/$line.2.name.tab
# Cat number of reads in raw files
  cat filtered_reads/$line.1.name.fastq | echo $((`wc -l`/4)) >> number_of_raw_reads.txt
# Removal of universal Illumina adapter and quality (q30) and length (45nt) trimming using Cutadapt
  cutadapt -b AGATCGGAAGAG -e 0.1 -j 8 -o filtered_reads/$line.1.name.adaptortrimmed.fastq -p filtered_reads/$line.2.name.adaptortrimmed.fastq filtered_reads/$line.1.name.fastq filtered_reads/$line.2.name.fastq # Removal of Universal Illumina adapter, at most 10% errors, using 8 cores
    rm filtered_reads/$line.1.name.fastq
    rm filtered_reads/$line.2.name.fastq
  cutadapt -q 30 -m 45 -j 8 -o filtered_reads/$line.1.q30.l45.fastq -p filtered_reads/$line.2.q30.l45.fastq filtered_reads/$line.1.name.adaptortrimmed.fastq filtered_reads/$line.2.name.adaptortrimmed.fastq # Filtering quality (-q) and length (-m), using 8 cores
    rm filtered_reads/$line.1.name.adaptortrimmed.fastq
    rm filtered_reads/$line.2.name.adaptortrimmed.fastq
# Cat number of reads in cleaned files
  cat filtered_reads/$line.1.q30.l45.fastq | echo $((`wc -l`/4)) >> number_of_cleaned_reads.txt
# Compress files
  pigz -p 24 $line.1.fastq
  pigz -p 24 $line.2.fastq
  pigz -p 24 filtered_reads/$line.1.q30.l45.fastq
  pigz -p 24 filtered_reads/$line.2.q30.l45.fastq
# Group file names, number of reads in raw files and number of reads in cleaned files
  paste $SAMPLELIST number_of_raw_reads.txt number_of_cleaned_reads.txt | column -s $'\t' -t > stats.tab
done

rm number_of_raw_reads.txt
rm number_of_cleaned_reads.txt

  set -o errexit   # abort on nonzero exitstatus
  set -o nounset   # abort on unbound variable
  set -o pipefail  # don't hide errors within pipes
#set -o xtrace          # Trace the execution of the script (debug)
