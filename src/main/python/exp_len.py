#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 3 :
    filename = sys.argv[1]
    refFlat_filename = sys.argv[2]
else:
    print("usage: python exp_len.py refSeq_MLE_output.tab known.gpd")
    print("or ./exp_len.py refSeq_MLE_output.tab known.gpd")
    sys.exit(1)
################################################################################

file = open(filename,'r')
row_i = 3 
dt = {}
for line in file:
    ls=line.strip().split('\t')
    dt[ ls[0] ] = ls[1] 
file.close()

################################################################################
used_set=set()
ref=open(refFlat_filename,'r')
len_dt={}
for refline in ref:
    refline_list=refline.strip().split()

    exon_start_list=refline_list[9].strip(',').split(',')
    exon_end_list=refline_list[10].strip(',').split(',')

    L = 0
    i=0
    for start in exon_start_list:
        start =int(start)
        end = int(exon_end_list[i])
        L += (end - start)
        i += 1
    if refline_list[1] in used_set:
        continue
    else:
        used_set.add(refline_list[1])
    if dt.has_key(refline_list[1]):
        print refline_list[0] + "\t" + refline_list[1] + "\t" + str(L) + "\t" + str(dt[refline_list[1]])
    else:
        print refline_list[0] + "\t" + refline_list[1] + "\t" + str(L) + "\t" + "0"

################################################################################
