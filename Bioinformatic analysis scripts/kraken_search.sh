#! /bin/bash
#$ -cwd
#$ -V
#$ -S /bin/bash

#Taxonomic attribution using Kraken2 and Bracken
#
#Author: Sarah François
#
#Date of last update: 2019-07-30
#
#USAGE: ./taxonomic-attribution-kraken SAMPLE-LIST
#
#EXAMPLE: ./taxonomic-attribution-kraken sample_list.txt

# SET VARIABLES
SAMPLELIST=$1
DBNAME=RDP

# DATA FORMATTING
  # fastq to fasta conversion
  cat $SAMPLELIST | while read line
  do
    seqtk seq -a /filtered_reads/$line.1.q30.l45.fastq.gz > $line.1.q30.fasta
    seqtk seq -a /filtered_reads/$line.2.q30.l45.fastq.gz > $line.2.q30.fasta
    cat $line.1.q30.fasta $line.2.q30.fasta > $line.q30.fasta
    pigz -p 24 $line.q30.fasta
  done
rm *.1.q30.fasta
rm *.2.q30.fasta

# CREATE DIRECTORY FOR RESULTS
mkdir Kraken-results

# KRAKEN-CLASSIFY A SAMPLE AND GENERATE REPORT FILES
cat $SAMPLELIST | while read line
do
  # TAXONOMIC ASSIGNMENT
      kraken2 --db $DBNAME --threads 24 --use-mpa-style --report Kraken-results/$line.mpa.report.kraken.tab --gzip-compressed $line.q30.fasta.gz > Kraken-results/$line.mpa.kraken.tab
      kraken2 --db $DBNAME --threads 24 --report Kraken-results/$line.report.kraken.tab --gzip-compressed $line.q30.fasta.gz > Kraken-results/$line.kraken.tab
  # BRACKEN [ABUNDANCE ESTIMATION] # -t = minimal number of reads per taxon threshold
      bracken -d $DBNAME -i Kraken-results/$line.report.kraken.tab -o Kraken-results/$line.brack.G.tab -r 150 -l 'G' -t 1 # classification at the genus scale
      sort -k 2n Kraken-results/$line.brack.G.tab > Kraken-results/$line.bracken.G.tab
      rm Kraken-results/$line.brack.G.tab
      bracken -d $DBNAME -i Kraken-results/$line.report.kraken.tab -o Kraken-results/$line.brack.F.tab -r 150 -l 'F' -t 1 # classification at the family scale
      sort -k 2n Kraken-results/$line.brack.F.tab > Kraken-results/$line.bracken.F.tab
      rm Kraken-results/$line.brack.F.tab
  # Files compression
      pigz -p 24 Kraken-results/$line.mpa.kraken.tab
      pigz -p 24 Kraken-results/$line.kraken.tab
      pigz -p 24 Kraken-results/$line.mpa.report.kraken.tab
      pigz -p 24 Kraken-results/$line.report.kraken.tab
done

# Move files into an appropriate folder
mkdir Kraken-results/mpa-report
mv Kraken-results/*.mpa.kraken.tab.gz Kraken-results/mpa-report/
mkdir Kraken-results/mpa-report-kraken
mv Kraken-results/*.mpa.report.kraken.tab.gz Kraken-results/mpa-report-kraken/
mkdir Kraken-results/report-kraken
mv Kraken-results/*.report.kraken.tab.gz Kraken-results/report-kraken/
mkdir Kraken-results/kraken
mv Kraken-results/*.kraken.tab.gz Kraken-results/kraken/
mkdir Kraken-results/report-kraken-bracken
mv Kraken-results/*.report.kraken_bracken.tab Kraken-results/report-kraken-bracken/

  set -o errexit   # abort on nonzero exitstatus
  set -o nounset   # abort on unbound variable
  set -o pipefail  # don't hide errors within pipes
#set -o xtrace          # Trace the execution of the script (debug)
