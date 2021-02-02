#! /usr/bin/env bash
#
# Author: Sarah François
# Date of last update: 2019-10-07
#
# Build KRAKEN2 + BRACKEN RDP database
#
# Usage: ./kraken_builddb.sh
#
# EXAMPLE: ./kraken_builddb.sh
#

# BUILD KRAKEN2 RDP DATABASE
kraken2-build --special rdp --db RDP --threads 24

# Build Bracken DATABASE
bracken-build -x kraken2/ -d RDP -t 24 -k 35 -l 150 -t 24
  # -k = kmer length
  # -l read length

set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes
#set -o xtrace          # Trace the execution of the script (debug)
