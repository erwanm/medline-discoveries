#!/usr/bin/env python3

#import sys
from os import listdir, mkdir
from os.path import isfile, join, isdir
from collections import defaultdict 
#import copy
#import spacy
#import scispacy


import sys, getopt


PROG_NAME = "extract-discoveries.py"

SEPARATOR_CONCEPT_ID_TYPE = "@"


COL_YEAR = 1
COLS_CONCEPT_INDIV = [ 2 ]
COL_COUNT_INDIV = 3 
COLS_CONCEPT_JOINT = [ 2, 3 ]
COL_COUNT_JOINT = 4 

OUTPUT_DETAILS = False
MIN_RATIO = None

def usage(out):
    print("Usage: "+PROG_NAME+" [options] <input file> <start year> <end year> <window size> <output file>",file=out)
    print("",file=out)
    print("  Reads the 'full' count file <input file> between <start year> and <end year> (inclusive) then",file=out)
    print("  calculates for every concept or pair of concepts the year where the increase between previous",file=out)
    print("  and next frequency is the highest. Frequency is taken as the sum over <window size> years.",file=out)
    print("  Notes:",file=out)
    print("    - The first <window_size> years from <start year> are ignored for the max ratio.",file=out)
    print("    - The last <window_size> years before <end year> are not ignored for the max ratio, but ",file=out)
    print("      they are at a disadvantage since the 'future' frequency is truncated.",file=out)
    print("",file=out)
    print("  Options")
    print("    -h: print this help message.",file=out)
    print("    -j: joint frequency  instead of individual frequency",file=out)
    print("    -p: PTC <concept>@<type> format as input; ignore the type.",file=out)
    print("    -d: details, add two column containing the sequence of raw fraquency and cumulated frequency.",file=out)
    print("    -m <min ratio>: keeps only concepts (or pairs) with a ratio at least equal to <min ratio> ",file=out)
    print("",file=out)


def convert_year_map_to_array(m, start, end):
    a = [0] * (end-start+1)
    for y,f in m.items():
        a[y-start] = f
    return a

# main

ptc_format = False
col_year = COL_YEAR -1
cols_concept = [ v-1 for v in COLS_CONCEPT_INDIV ]
col_freq = COL_COUNT_INDIV -1
try:
    opts, args = getopt.getopt(sys.argv[1:],'hjpdm:')
except getopt.GetoptError:
    usage(sys.stderr)
    sys.exit(2)
for opt, arg in opts:
    if opt == "-h":
        usage(sys.stdout)
        sys.exit()
    elif opt == "-j":
        cols_concept = [ v-1 for v in COLS_CONCEPT_JOINT ]
        col_freq = COL_COUNT_JOINT -1
    elif opt == '-p':
        ptc_format = True
    elif opt == '-d':
        OUTPUT_DETAILS = True
    elif opt == '-f':
        MIN_RATIO = float(arg)



#print("debug args after options: ",args)

if len(args) != 5:
    usage(sys.stderr)
    sys.exit(2)

input_file = args[0]
start_year = int(args[1])
end_year = int(args[2])
window_size = int(args[3])
output_file = args[4]

# count_map[concept][year] = frequency
count_map = defaultdict(lambda: defaultdict(int))
with open(input_file) as infile:
    for line in infile:
        cols = line.rstrip().split("\t")
        year = int(cols[col_year])
        if year >= start_year and year <= end_year:
            freq = int(cols[col_freq])
            concepts = []
            if ptc_format:
                for i in cols_concept:
                    val = cols[i]
                    sep_pos = val.rfind(SEPARATOR_CONCEPT_ID_TYPE)
                    if sep_pos != -1:
                        #                ptc_category = concept[sep_pos+1:]
                        val = val[0:sep_pos]
                    else:
                        print("Warning: PTC format separator not found in '"+val+"'",file=sys.stderr)
                        concepts.append(val)
                        concept = '\t'.join(concepts)
            else:
                concept = '\t'.join([ cols[i] for i in cols_concept ])
            count_map[concept][year] += freq

with open(output_file,"w") as outfile:
    for concept,map_by_year in count_map.items():
#        print(map_by_year)
        d = convert_year_map_to_array(map_by_year, start_year, end_year)
#        print(d)
        # cumulated[i] = [ past, future, ratio ]
        cumulated = []
        max_ratio_year_index = None
        total = 0
        past = 0
        future = 0
        # init future freq with first N years:
        for i in range(0, window_size):
            future += d[i]
        for i in range(len(d)):
            total += d[i]
            if i > 0:
                past += d[i-1]
                future -= d[i-1]
                if i+window_size-1 < end_year-start_year+1: # not done for i=0 since already done during init
                    future += d[i+window_size-1]
                if i > window_size:
                    past -= d[i-window_size-1]
            ratio = None
            if past > 0:
                ratio = future / past
            cumulated.append([past, future, ratio ])
            if i >= window_size and ratio is not None and (max_ratio_year_index is None or ratio > cumulated[max_ratio_year_index][2]):
                max_ratio_year_index = i 
        ratio = None
        if max_ratio_year_index is not None:
            ratio = cumulated[max_ratio_year_index][2]
        if MIN_RATIO is None or (ratio is not None and ratio >= MIN_RATIO):
            if max_ratio_year_index is None:
                str_y = "NA"
            else:
                str_y = str(start_year+max_ratio_year_index)
            if ratio is None:
                str_ratio = "NA"
            else:
                str_ratio = str(ratio)
            outfile.write("%s\t%s\t%s\t%d" % (concept, str_y, str_ratio, total))
            if OUTPUT_DETAILS:
                raw_freq_str = ' '.join([str(v) for v in d])
                cumulated_freq_str = ' '.join([ str(c[0])+','+str(c[1]) for c in  cumulated ])
                outfile.write("\t%s\t%s" % (raw_freq_str, cumulated_freq_str))
            outfile.write("\n")
