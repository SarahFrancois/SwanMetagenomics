#! /usr/bin/env bash
#
# Author: Sarah François, CC BY-NC-SA
# Date of last update: 2019-06-12
#
# Usage: SCRIPTNAME [OPTIONS]... [ARGUMENTS]...
#
#
# OPTIONS
#	-h Print the help message
#
# EXAMPLE
# 

#Set magic variables for current file, basename, and directory
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
__file="${__dir}/$(basename "${BASH_SOURCE[0]}")"
__base="$(basename ${__file} .sh)"
__root="$(cd "$(dirname "${__dir}")" && pwd)" # <-- change this as it depends on your app
arg1="${1:-}"

# Help 
if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` [somestuff]"
  exit 0
fi



  set -o errexit   # abort on nonzero exitstatus
  set -o nounset   # abort on unbound variable
  set -o pipefail  # don't hide errors within pipes
#set -o xtrace          # Trace the execution of the script (debug)

