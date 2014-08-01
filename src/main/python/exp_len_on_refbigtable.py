#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 3 :
    ref_filename = sys.argv[1]
    tag_filename = sys.argv[2]
else:
    print("usage: H1_exp/multiexon_refFlat.txt_positive_exp_len known_intact_SM.fa.bestpsl.gpd_refFlat.txt")
    print("or ")
    sys.exit(1)
################################################################################
tag_s= set()
tag = open(tag_filename,'r')
for line in tag:
    tag_s.add(line.strip().split("\t")[1])
tag.close()

ref = open(ref_filename,'r')
for line in ref:
    ls = line.strip().split("\t")
    i = "0"
    if ls[1] in tag_s:
        i = "1"
    ls.append(i)
    print "\t".join(ls)
ref.close()
################################################################################
    
