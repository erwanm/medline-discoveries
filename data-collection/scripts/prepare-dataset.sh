#!/bin/bash

minFreq=100
optPTC=""

function usage {
    echo
    echo "Usage: $0 <indiv freq dir> <joint freq dir> <output dir>" 
    echo 
    echo "  - Requires tdc-tools/code in PATH."
    echo "  - the output dir can be the one which contains the 2 input dirs (subdirs created)" 
    echo
    echo "  Options:"
    echo "    -h: print this help message."
    echo "    -p: PTC data with <concept>@<type> format"
    echo "    -m <min freq> default: 100"
    echo ""
}





while getopts 'hpm:' option ; do
    case $option in
        "h" ) usage
              exit 0;;
	"p" ) optPTC="-p";;
	"m" ) minFreq="$OPTARG";;
        "?" )
            echo "Error, unknow option." 1>&2
            usage 1>&2
                exit 1
    esac
done
shift $(($OPTIND - 1)) # skip options already processed above                                                                                                            
if [ $# -ne 3 ]; then
    echo "Error: 3 arguments expected." 1>&2
    echo 1>&2
    usage 1>&2
    exit 1
fi


indivFreqDir="$1"
jointFreqDir="$2"
outputDir="$3"

rm -f "$indivFreqDir"/*.err "$jointFreqDir"/*.err

acrossDir="$outputDir/across-all-years"
[ -d "$acrossDir" ] || mkdir "$acrossDir"

# generate indiv data across years with filtered min freq version
if [ -z "$optPTC" ]; then
    sum-freq-over-years.py "$indivFreqDir" 0 3000 "$acrossDir"/indiv.full
else
    tmpPTC=$(mktemp --tmpdir=/tmp prepare-dataset.XXXXXXXX)
    sum-freq-over-years.py "$indivFreqDir" 0 3000 "$tmpPTC"
    ptc-aggregate-across-types.py "$tmpPTC" "$acrossDir"/indiv.full
    cat "$tmpPTC.total" >"$acrossDir"/indiv.full.total
    rm -f "$tmpPTC"
fi
echo $minFreq | filter-column.py -m "$acrossDir"/indiv.full 2 >"$acrossDir"/indiv.min${minFreq}

# generate indiv data by year with filtered min freq version
tmpCat=$(mktemp --tmpdir=/tmp prepare-dataset.XXXXXXXX)
cat "$indivFreqDir"/???? >"$tmpCat"
cat "$acrossDir"/indiv.min${minFreq} | cut -f 1 | filter-column.py $optPTC "$tmpCat" 2 >"$outputDir"/indiv.full.min${minFreq}
rm -f $tmpCat
cat "$indivFreqDir"/*.total >"$outputDir"/indiv.full.total

# generate joint data across years with filtered min freq version
if [ -z "$optPTC" ]; then
    sum-freq-over-years.py -j "$jointFreqDir" 0 3000 "$acrossDir"/joint.full
else
    tmpPTC=$(mktemp --tmpdir=/tmp prepare-dataset.XXXXXXXX)
    sum-freq-over-years.py -j "$jointFreqDir" 0 3000 "$tmpPTC"
    ptc-aggregate-across-types.py -j "$tmpPTC" "$acrossDir"/joint.full
    rm -f "$tmpPTC"
fi
echo $minFreq | filter-column.py -m "$acrossDir"/joint.full 3 >"$acrossDir"/joint.min${minFreq}

# generate joint data by year with filtered min freq version
tmpCat=$(mktemp --tmpdir=/tmp prepare-dataset.XXXXXXXX)
cat "$jointFreqDir"/???? >"$tmpCat"
cat "$acrossDir"/joint.min${minFreq} | cut -f 1,2 | filter-column.py  $optPTC $tmpCat 2,3 >"$outputDir"/joint.full.min${minFreq}
rm -f $tmpCat

