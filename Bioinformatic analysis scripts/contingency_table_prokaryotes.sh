#! /usr/bin/env bash
#
# Author: Sarah François
# Date of last update: 2019-10-07
#
# Usage: ./contingency_table_prokaryotes.sh SAMPLE_LIST BRACKEN_FILE_NAME > contingency_table_prokaryotes.tab
#
# EXAMPLE: ./contingency_table_prokaryotes.sh sample_list.txt bracken.F.tab> contingency_table_prokaryotes.tab
#

# SET VARIABLES
SAMPLELIST=$1

echo -n "Taxonomy" ; for i in `cat $SAMPLELIST` ; do echo -n "@"$i"" ; done | tr '@' '\t'

echo ""

for i in `cat $SAMPLELIST` ; do cat "$i"."$2" ; done  | awk '{print $1}' | sort | uniq > Taxonomy.tmp

for i in `cat $SAMPLELIST` ; do for j in `cat Taxonomy.tmp` ; do if [ `grep -w "$j" "$i"."$2" | awk '{print $1}' | uniq` ] ; then grep -w "$j" "$i"."$2" | awk '{print sum+=$2}' | tail -1 ; else echo "0" ; fi ; done > "$i".tmp ; done

list=`ls *.tmp | grep -f $SAMPLELIST`

nb=`wc -l $SAMPLELIST | awk '{print $1+1}'`

pr -tmw5000 Taxonomy.tmp `echo "$list"` | awk '{for (i=1;i<='"$nb"';i++) {printf $i"@"} {print ""}}' | sed -e 's/@$//g' | tr '@' '\t'

for i in `cat $SAMPLELIST` ; do rm "$i".tmp ; done ; rm Taxonomy.tmp

  set -o errexit   # abort on nonzero exitstatus
  set -o nounset   # abort on unbound variable
  set -o pipefail  # don't hide errors within pipes
#set -o xtrace          # Trace the execution of the script (debug)
