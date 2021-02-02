#! /usr/bin/env bash
#
# Author: Sarah François
# Date of last update: 2019-10-07
#
# Usage: ./mapping.sh CONTIGS SAMPLELIST
#
# EXAMPLE: ./mapping.sh SPAdes-assembly/contigs.fasta sample_list.txt
#

# SET VARIABLES
CONTIGS=$1
SAMPLELIST=$2

# CREATE DIRECTORY FOR RESULTS
mkdir Mapping

# DATABASE CREATION USING BOWTIE2
bowtie2-build -f --threads 4 $CONTIGS Mapping/VIRUS.DB

# MAPPING OF PAIRED END READS (BOWTIE2) + CALCUL OF NUMBER OF READS ASSIGNED TO EACH CONTIG (bbmap)
cat $SAMPLELIST | while read line
do
  bowtie2 --threads 24 -x VIRUS.DB -1 /filtered_reads/$line.1.q30.l45.fastq.gz -2 /filtered_reads/$line.2.q30.l45.fastq.gz -S Mapping/$line.virus.sam --no-unal --end-to-end --al-conc Mapping/$line.mapped.virus.fastq
  bbmap/pileup.sh in=Mapping/$line.virus.sam out=Mapping/$line.virus.cov.tab
done

# FILE GENERATION FOR CONTINGENCY TABLE CREATION
cat Mapping/*.virus.sam > Mapping/mapped_reads.sam

  set -o errexit   # abort on nonzero exitstatus
  set -o nounset   # abort on unbound variable
  set -o pipefail  # don't hide errors within pipes
#set -o xtrace          # Trace the execution of the script (debug)
