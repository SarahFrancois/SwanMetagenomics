#! /usr/bin/env bash
#
# Author: Sarah François
# Date of last update: 2019-10-07
#
# Usage: ./de_novo_assembly.sh SAMPLE_LIST
#
# EXAMPLE: ./de_novo_assembly.sh sample_list.txt
#

# SET VARIABLES
SAMPLELIST=$1

# CREATE DIRECTORY FOR RESULTS
mkdir SPAdes-assembly

# DE NOVO ASSEMBLY USING SPAdes
cat $SAMPLELIST | while read line
do
  spades.py -k 21,33,55,77 -t 64 --careful -1 /filtered_reads/$line.1.q30.l45.fastq.gz -2 /filtered_reads/$line.2.q30.l45.fastq.gz -o SPAdes-assembly/$line.Assembly
done

  set -o errexit   # abort on nonzero exitstatus
  set -o nounset   # abort on unbound variable
  set -o pipefail  # don't hide errors within pipes
#set -o xtrace          # Trace the execution of the script (debug)
