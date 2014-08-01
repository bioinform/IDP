#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 3:
    psl_filename = sys.argv[1]
    skip_Nine = int(sys.argv[2])
else:
    print("usage:hist_psl.py psl_file skip_Nine")
    print("or ")
    sys.exit(1)
################################################################################
def process_temp_list(temp_list):
    ref_stat = 0
    ref_ls = []
    for result_ls in temp_list:
        stat = float(result_ls[0])/float(result_ls[10])
        if stat > ref_stat:
            ref_stat = stat
            ref_ls = result_ls
    return ref_ls
################################################################################
psl = open(psl_filename,'r')
temp_list = []
Qname=""
i=0
for line in psl:
    if i<skip_Nine:
        i+=1
        continue
    ls = line.strip().split('\t')
    if Qname == ls[9]:
        temp_list.append(ls)
    else:
        if not Qname =="":
            result_ls = process_temp_list(temp_list)
            print '\t'.join(result_ls)
        Qname = ls[9]
        temp_list = [ls]
    i+=1
result_ls = process_temp_list(temp_list)
print '\t'.join(result_ls)
psl.close()
################################################################################

