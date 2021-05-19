#!/bin/bash

set -e 

if [ $# -ne 3 ]; then
   echo "need 3 input arguments, (1) input file, (2) superpopulations (e.g. 'AMR EUR ...') and (3) output prefix"
   exit 1
fi

# checks if the supplied superpopulations match those in the file
head -n 1 $1 | awk -v sp="$2" 'BEGIN{split(sp,array)}{{for(i=3;i<=NF;i++){if($i =! array[i-2]){exit 1}}}}'

# writes keep files for the different ancestries with stringent cutoff
awk -v sp="$2" -v out="$3" 'BEGIN{split(sp,array)}{if(NR>1){{for(i=3;i<=NF;i++){if($i > 0.995){print $1"\t"$2 > out "." array[i-2]".keep"}}}}}' $1