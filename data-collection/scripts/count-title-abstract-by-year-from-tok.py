#!/usr/bin/env python3

#import sys
from os import listdir, mkdir
from os.path import isfile, join, isdir
from collections import defaultdict 
#import copy
#import spacy
#import scispacy


import sys, getopt


PROG_NAME = "count-title-abstract-by-year-from-tok.py"


def usage(out):
    print("Usage: "+PROG_NAME+" [options] <input file> <output file>",file=out)
    print("",file=out)
    print("  Reads a tokenized .tok file (TDC data), counts for every document the number of tokens in",file=out)
    print("  the title, abstract and in total.",file=out)
    print("  Options")
    print("    -h: print this help message.",file=out)
    print("",file=out)



# main

try:
    opts, args = getopt.getopt(sys.argv[1:],'h')
except getopt.GetoptError:
    usage(sys.stderr)
    sys.exit(2)
for opt, arg in opts:
    if opt == "-h":
        usage(sys.stdout)
        sys.exit()



#print("debug args after options: ",args)

if len(args) != 2:
    usage(sys.stderr)
    sys.exit(2)

input_file = args[0]
output_file = args[1]

# count_map[pmid][part] = token count, where 'part' is abstract or title 
count_tokens = defaultdict(lambda: defaultdict(int))
count_sent = defaultdict(int)
year = None
with open(input_file) as infile:
    for line in infile:
#        print(line)
        cols = line.rstrip("\n").split("\t")
        pmid = cols[0]
        this_year = cols[1]
        if year is None:
            year = this_year
        else:
            if this_year != year:
                print("Warning: inconsistent year found in ",input_file,":",this_year," (expected ",year,")",file=sys.stderr)
        part = cols[2]
        tokens = cols[5]
        count = len(tokens.split(" "))
        count_tokens[pmid][part] += count
        count_sent[pmid] += 1

with open(output_file,"w") as outfile:
    for pmid,m in count_tokens.items():
        sent = count_sent[pmid]
        count_title = m.get('title',0)
        count_abstract = m.get('abstract',0)
        outfile.write("%s\t%s\t%d\t%d\t%d\t%d\n" % (year, pmid, sent, count_title, count_abstract,count_title + count_abstract))


