#! /usr/bin/env bash
#
# Author: Sarah François
# Date of last update: 2019-10-07
#
# Usage: ./diamond_search.sh CONTIGS
#
# EXAMPLE: ./diamond_search.sh SPAdes-assembly/contigs.fasta
#

# SET VARIABLES
CONTIGS=$1
i=${1%.gz}
i=${i%.fnas}
i=${i%.fastq}
i=${i%.fasta}
j=`echo $i | sed -e 's/.*\///g'`

# Diamond blast against viral database
diamond blastx -q $CONTIGS -d VIRUSES/VIRUSES.dmnd -o "$j".tmp -f 6 qseqid full_qseq -k 1 -c 1 -b 8 -e 0.0001 -v # diamond blast
awk  -F '\t' '{print ">"$1"\n"$2}' "$j".tmp > "$j".tmp2 # diamond positive hits query sequences, fasta format
# seqtk subseq "$1" "$j".tmp > "$j".tmp2 # diamond positive hits query sequences, fastq format

# Diamond blast against nr database
diamond blastx -q "$j".tmp2 -d NR/NR.dmnd -o "$j".tmp3 -f 6 qseqid full_qseq qlen stitle pident evalue sseqid length mismatch gapopen qstart qend sstart send bitscore staxids -k 1 -c 1 -b 8 -e 0.0001 -v

# Formatting and creation of diamond blast results files
  # Formatting files to print taxids from diamond blast output list
cat "$j".tmp3 | sort -t$'\t' -k16,16n > "$j".tmp4 ; awk -F"\t" '{print $16}' "$j".tmp4 > "$j".tmp5 # taxids list
cat "$j".tmp5 | /usr/local/bin/ete3 ncbiquery --info | grep -v "^#" > "$j".tmp6 # print taxonomy from taxids
  # Add taxonomy to diamond blastx results file
join -t$'\t' -a1 -e'NA' -116 -21 "$j".tmp4 "$j".tmp6 -o 1.1,1.2,1.3,1.4,2.2,2.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16 > "$j".tmp7
  # Creation of diamond blastx results in tabular format, and their corresponding sequences in fasta format
cat "$j".tmp7 | grep "Viruses,\|virus" | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18}' > "$j".tmp8
awk -F"\t" '$6!="NA"' $j.tmp8 > $j.tmp9
awk -F"\t" '{print $1"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18}' $j.tmp9 > $j.tmp10
  head=( qseqid qlen sseqid sseqidlong stitle taxonomy pident evalue length mismatch gapopen qstart qend sstart send bitscore staxids )
  ( IFS=$'\t'; echo "${head[*]}"; cat $j.tmp10 ) > $j.blastx.tab # tabular file with headers
awk -F"\t" '{print ">"$1"\n"$2}' $j.tmp9 > $j.blastx.fasta # fasta file

rm "$j".tmp*

  set -o errexit   # abort on nonzero exitstatus
  set -o nounset   # abort on unbound variable
  set -o pipefail  # don't hide errors within pipes
#set -o xtrace          # Trace the execution of the script (debug)
