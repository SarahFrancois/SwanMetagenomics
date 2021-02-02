#! /bin/bash
#$ -cwd
#$ -V
#$ -S /bin/bash

# Author: Sarah Francois

# Impact of host age on microbial communities in the wild: insights from a waterbird population - pipeline

# Custom scripts used for data analysis of metagenomic datasets
  # Objectives: Viral genomes reconstruction and contingency tables creation
  # Input files: (1) List of Illumina HiSeq 4000 sequencing paired end reads in fastq format and (2) list of sequenced samples (one sample per row)
  # Example of input files: (1) UTMF0973.R1.fastq.gz UTMF0973.R2.fastq.gz UTMF0975.R1.fastq.gz UTMF0975.R2.fastq.gz (2) UTMF0973 UTMF0975

# Quality filtering
./quality_filtering.sh sample_list.txt
  # filters high-throughput sequencing reads >q30 quality and removes reads >45 bases length using Cutadapt

# De novo assembly
./de_novo_assembly.sh sample_list.txt
  # de novo assembly of cleaned reads using SPAdes

# Taxonomic assignment
  # Viruses - DIAMOND
    # Build DIAMOND databases
    ./diamond_builddb.sh
    # Taxonomic assignment using DIAMOND
    ./diamond_search.sh SPAdes-assembly/contigs.fasta
  # Prokaryotes - KRAKEN2 + BRACKEN
    # Build KRAKEN RDP database
    ./kraken_builddb.sh
    # Taxonomic assignment using KRAKEN2 and refinement using BRACKEN
    ./kraken_search.sh

# Mapping of viral contigs
./mapping.sh SPAdes-assembly/contigs.fasta sample_list.txt

# Contingency table creation for viral contigs
./contingency_table_viruses.sh sample_list.txt # viruses
./contingency_table_prokaryotes.sh sample_list.txt # prokaryotes
