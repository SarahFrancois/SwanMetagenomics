#! /usr/bin/env bash
#
# Author: Sarah François
# Date of last update: 2019-10-07
#
# Build DIAMOND NCBI nr and refseq virus databases
#
# Usage: ./diamond_builddb.sh
#
# EXAMPLE: ./diamond_builddb.sh
#

# Requirements for DIAMOND databases creation: fasta datasets, NCBI protein accession numbers and NCBI taxonomy
  # Download of nr and refseq viral datasets
    wget ftp://ftp.ncbi.nlm.nih.gov/blast//db/FASTA/nr.gz # nr dataset
    wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/viral.1.protein.faa.gz # refseq viral dataset
  # Download and uncompress NCBI protein accession numbers to taxon ids
    mkdir TAXONOMY
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz -P TAXONOMY
    gunzip TAXONOMY/prot.accession2taxid.gz
  # Download and uncompress NCBI taxonomy
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip -P TAXONOMY
    unzip TAXONOMY/taxdmp.zip

# Build databases:
  # nr database
    mkdir NR
    diamond makedb --in nr.gz --db NR/NR --taxonmap TAXONOMY/prot.accession2taxid --taxonnodes TAXONOMY/nodes.dmp --taxonnames TAXONOMY/names.dmp --threads 48
  # refseq viral database
    mkdir VIRUSES
    diamond makedb --in viral.1.protein.faa.gz --db VIRUSES/VIRUSES --taxonmap TAXONOMY/prot.accession2taxid --taxonnodes TAXONOMY/nodes.dmp --taxonnames TAXONOMY/names.dmp --threads 48

# Removal of nr.gz and viral.1.protein.faa.gz files
  rm nr.gz
  rm viral.1.protein.faa.gz

  set -o errexit   # abort on nonzero exitstatus
  set -o nounset   # abort on unbound variable
  set -o pipefail  # don't hide errors within pipes
#set -o xtrace          # Trace the execution of the script (debug)
