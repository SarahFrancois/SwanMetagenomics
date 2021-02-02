# SwanMetagenomics

Author: Sarah Francois

Custom scripts used for data analysis of viral and 16S metagenomic datasets
  
## Bioinformatic pipeline
  Objectives: Viral genomes reconstruction and contingency tables creation
  Input files: (1) List of Illumina HiSeq 4000 sequencing paired end reads in fastq format and (2) list of sequenced samples (one sample per row)
  Example: (1) UTMF0973.R1.fastq.gz UTMF0973.R2.fastq.gz UTMF0975.R1.fastq.gz UTMF0975.R2.fastq.gz (2) UTMF0973 UTMF0975

### Required sofwares
- cutadapt
- SPAdes
- DIAMOND
- kraken2
- bracken
- BOWTIE2
- bbmap
- seqtk

### Scripts
 
 #### Quality filtering
- quality_filtering.sh
 
 #### De novo assembly
- de_novo_assembly.sh
 
 #### Taxonomic assignment
  ##### Viruses - DIAMOND
- diamond_builddb.sh  # Build databases
- diamond_search.sh SPAdes-assembly/contigs.fasta # Taxonomic assignment
  ##### Prokaryotes - KRAKEN2 + BRACKEN
- kraken_builddb.sh # Build database
- kraken_search.sh # Taxonomic assignment

#### Mapping of viral contigs
- mapping.sh

#### Contingency tables creation
- contingency_table_viruses.sh # viruses
- contingency_table_prokaryotes.sh # prokaryotes


## Statistical pipeline 
  Objectives: Diversity and differential prevalence and abundance analysis of prokaryotic and viral communities
  Input files: (1) Contingency table in tabular format and (2) corresponding samples metadata table in tabular format
  Example: (1) contingency_table.tab (2) metadata.tab
  
### Required sofwares
- R
- RStudio

### Script
swan-metagenomics.R
