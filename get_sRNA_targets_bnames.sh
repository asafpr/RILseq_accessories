#!/usr/bin/env sh
# Get the targets of a sRNA and write their bnumbers in a file
# $1 - sRNA name
# $2 - interactions file
# $3 - NC_000913.ptt.gz path
grep -w $1 $2 | cut -f 3 | tr "." "\n"  | sort | uniq > list_of_${1}_targets_names.txt
grep -w $1 $2 | cut -f 4 | tr "." "\n"  | sort | uniq >> list_of_${1}_targets_names.txt
sort list_of_${1}_targets_names.txt | uniq | grep -v $1 > list_of_${1}_targets_names.txt2
zcat $3 | grep -w -f list_of_${1}_targets_names.txt2 | cut -f 6 
