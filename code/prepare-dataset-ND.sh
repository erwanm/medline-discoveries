#!/bin/bash

#set -x

minFreq=100
optAddTerm=""
topPMI=100


function usage {
    echo
    echo "Usage: $0 <full prepared dir> <targets file> <UMLS dir> <UMLS semantic groups file>" 
    echo 
    echo "  - input dir = result of prepare-dataset.sh, contains indiv.full.min100 and joint.full.min100"
    echo "  - resulting files created in the same directory"
    echo "  - Requires tdc-tools/code and medline-discoveries/code in PATH."
    echo "  - <UMLS sem groups file> UMLS group hierarchy file which can be downloaded from:"
    echo "      https://lhncbc.nlm.nih.gov/semanticnetwork/download/SemGroups.txt" 

    echo
    echo "  Options:"
    echo "    -h: print this help message."
    echo "    -p: PTC data: MeSH decriptors with prefix (for add-terms)"
    echo "    -m: MED data: MeSH descriptors (for add-terms)"
    echo "    -f <min freq> default: 100. input files named <indiv|joint>.full.min<min freq>"
    echo "    -t <n> number of top PMI/NPMI to extract from the ND 'across years' data."
    echo ""
}





while getopts 'hpmf:t:' option ; do
    case $option in
        "h" ) usage
              exit 0;;
	"p" ) optAddTerm="-M";;
	"m" ) optAddTerm="-m";;
	"f" ) minFreq="$OPTARG";;
	"t" ) topPMI="$OPTARG";;
        "?" )
            echo "Error, unknow option." 1>&2
            usage 1>&2
                exit 1
    esac
done
shift $(($OPTIND - 1)) # skip options already processed above                                                                                                            
if [ $# -ne 4 ]; then
    echo "Error: 4 arguments expected." 1>&2
    echo 1>&2
    usage 1>&2
    exit 1
fi


inputDir="$1"
targetsFile="$2"
umlsDir="$3"
semGroups="$4"

indivFullFile="$inputDir/indiv.full.min${minFreq}"
jointFullFile="$inputDir/joint.full.min${minFreq}"
indivNDFile="$inputDir/indiv.ND.min${minFreq}"
jointNDFile="$inputDir/joint.ND.min${minFreq}"
totalFile="$inputDir/full.total"

noYearDir="$inputDir/across-all-years"
noYearIndiv="$noYearDir/indiv.min${minFreq}"
noYearJoint="$noYearDir/joint.min${minFreq}"

# 1. filter ND cooccurrences in the "across years" data

#  1. METHOD A
tmpIndiv=$(mktemp --tmpdir=/tmp "$(basename $0).indiv.XXXXXXXXXX")
#    (a) collect cooccurrences with ND targets
cat "$targetsFile" | filter-column.py -u "$noYearJoint" 1,2 | cut -f 1,2 | tr '\t' '\n' | sort -u >$tmpIndiv
#    (b) collect co-occuring concepts (lincludes targets themselves)
cat  $tmpIndiv | filter-column.py "$noYearIndiv" 1 >"$noYearIndiv.ND"
#    (c) filter joint data with any pair of co-occuring concept
cat  $tmpIndiv | filter-column.py -i "$noYearJoint" 1,2 >"$noYearJoint.ND"
rm -f "$tmpIndiv"

#  1. METHOD B: using 'GROUP:ND' - ABANDONED, see my notes
#outputGroup="$noYearDir/group-ND.concepts"
#tmpJoint=$(mktemp --tmpdir=/tmp "$(basename $0).indiv.XXXXXXXXXX")
#cat "$noYearJoint" | grep "GROUP:ND"  > "$tmpJoint"
#calculate-association-measure.py -m pmi,npmi "$tmpJoint" "$noYearIndiv" "$noYearDir/indiv.full.total" "$outputGroup.pmi"
#rm -f "$tmpJoint"
##   filter positive PMI pairs
#echo 0 | filter-column.py -m "$outputGroup.pmi" 4 >"$outputGroup.pos-pmi"
#tmpIndiv=$(mktemp --tmpdir=/tmp "$(basename $0).indiv.XXXXXXXXXX")
#cut -f 1,2 "$outputGroup.pos-pmi" | tr '\t' '\n' | grep -v "GROUP:ND" | sort -u >$tmpIndiv
#cat  $tmpIndiv | filter-column.py "$noYearIndiv" 1 >"$noYearIndiv.ND.group-version"
#cat  $tmpIndiv | filter-column.py -i "$noYearJoint" 1,2 >"$noYearJoint.ND.group-version"
#rm -f "$tmpIndiv"

# 2. calculate PMI in the "across years" data
calculate-association-measure.py -m pmi,npmi "$noYearJoint.ND" "$noYearIndiv.ND" "$noYearDir/indiv.full.total" "$noYearJoint.ND.pmi"
#   filter positive PMI pairs
echo 0 | filter-column.py -m "$noYearJoint.ND.pmi" 4 >"$noYearJoint.ND.pos-pmi"
#   extract top N PMI and NPMI
sort -r -g -k4,4 "$noYearJoint.ND.pmi" | head -n $topPMI >"$noYearJoint.ND.pmi.top${topPMI}"
sort -r -g -k5,5 "$noYearJoint.ND.pmi" | head -n $topPMI >"$noYearJoint.ND.npmi.top${topPMI}"
# merge the two top lists
tmpBoth=$(mktemp --tmpdir=/tmp "$(basename $0).both.XXXXXXXXXX")
cat "$noYearJoint.ND.pmi.top${topPMI}" "$noYearJoint.ND.npmi.top${topPMI}" | sort -u >"$tmpBoth"
# add term and category for top N
ls "$tmpBoth" | add-term-from-umls.py $optAddTerm -i 1 -g "$semGroups" "$umlsDir" .term1
ls "$tmpBoth.term1" | add-term-from-umls.py $optAddTerm -i 2 -g "$semGroups" "$umlsDir" .term2
mv "$tmpBoth.term1.term2" "$noYearJoint.ND.top-pmi-npmi"
rm -f $tmpBoth $tmpBoth.term1

# 3. apply filter to "by year" data
cut -f 1 "$noYearIndiv.ND" | filter-column.py "$indivFullFile" 2 >"$indivNDFile"
cut -f 1,2 "$noYearJoint.ND" | filter-column.py "$jointFullFile" 2,3 >"$jointNDFile"
cut -f 1,2 "$noYearJoint.ND.pos-pmi" | filter-column.py "$jointFullFile" 2,3 >"$jointNDFile.pos-pmi"

