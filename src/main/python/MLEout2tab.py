#!/usr/bin/python

import sys
import os
if len(sys.argv) >= 2 :
    exp_filename = sys.argv[1]
else:
    print("usage: ")
    print("or ")
    sys.exit(1)
################################################################################
exp_dt = {}
exp_file = open(exp_filename,'r')
i=0
for line in exp_file:
    ls = line.strip().split('\t')
    if 1:
        if i%3 == 0:
            ls = line.strip().split('\t')
            gene_name = ls[0]
        elif i%3 == 1:
            ID_ls = line.strip().split('\t')
        elif i%3 == 2:
            exp_str_ls = line.strip().split('\t')
            j = 0
            for ID in ID_ls:
                if not exp_dt.has_key(gene_name):
                     exp_dt[gene_name] = {}
                if not exp_str_ls[j] == "NA":
                    exp_dt[gene_name][ID]= float(exp_str_ls[j])
                j += 1
        i+=1
exp_file.close()
#print exp_dt
################################################################################

for gene_name in exp_dt:
    for ID in exp_dt[gene_name]:
        if ID!="Total":
            print '\t'.join([ID, str(exp_dt[gene_name][ID]), str(exp_dt[gene_name]["Total"])] )
