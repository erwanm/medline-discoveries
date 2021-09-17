#!/bin/bash


ptc=""
if [ "$1" == "-p" ]; then
    ptc="-p"    
    shift
fi

if [ $# -ne 5 ]; then
    echo "Usage: $0 [-p] <targets file> <input indiv> <input joint> <output indiv> <output joint>" 1>&2
    exit 1
fi


targetsFile="$1"
inputIndiv="$2"
inputJoint="$3"
outputIndiv="$4"
outputJoint="$5"

tmpIndiv=$(mktemp --tmpdir=/tmp "$(basename $0).indiv.XXXXXXXXXX")
#echo "DEBUG: $tmpIndiv" 1>&2
cat "$targetsFile" | filter-column.py $ptc -u "$inputJoint" 2,3 | cut -f 2,3 | tr '\t' '\n' | sort -u >$tmpIndiv
cat  $tmpIndiv | filter-column.py $inputIndiv 2 >$outputIndiv
cat  $tmpIndiv | filter-column.py -i $inputJoint 2,3 >$outputJoint
