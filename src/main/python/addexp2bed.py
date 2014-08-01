#!/usr/bin/python

import sys
import os
col_str = "0,0,0"
index_ls = [0,1,2]
if len(sys.argv) >= 4 :
    exp_filename = sys.argv[1]
    gpd_filename = sys.argv[2]
    output_filename = sys.argv[3]
    if len(sys.argv)>=5:
        if sys.argv[4].upper() == "B":
            col_str = "128,0,255"
            index_ls = [1]
        elif sys.argv[4].upper() == "G":
            col_str = "0,255,0"
            index_ls = [0,2]
        elif sys.argv[4].upper() == "R":
            col_str = "255,0,128"
            index_ls = [1]
        elif sys.argv[4].upper() == "Y":
            col_str = "255,255,0"
            index_ls = [2]
else:
    print("usage: ")
    print("or ")
    sys.exit(1)
################################################################################
def design_col(score,range_start,range_end,col_str,index):
    ls = col_str.split(',')
    ls[index] = str( range_end - int((range_end - range_start) * float(score)/1000) )
    return ','.join(ls)
################################################################################
exp_dt = {}
exp_file = open(exp_filename,'r')
i=0
for line in exp_file:
    line=line.strip()
    if line == '':
        i = 0
        continue
    else:
        if i == 0:
            ls = line.strip().split('\t')
            gene_name = ls[0]
        elif i == 1:
            ID_ls = line.strip().split('\t')
            del ID_ls[-1]
        elif i == 2:
            exp_str_ls = line.strip().split('\t')
            del exp_str_ls[-1]
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
score_dt = {}
for gene_name in exp_dt:
    score_dt[gene_name] = {}
    total = 0
    for ID in exp_dt[gene_name]:
        total += exp_dt[gene_name][ID]
    if total == 0:
        continue
    for ID in exp_dt[gene_name]:
        score_dt[gene_name][ID] = int( exp_dt[gene_name][ID]*1000/total )

#print score_dt
################################################################################
gpd_file = open(gpd_filename,'r')
output = open(output_filename,'w')
for line in gpd_file:
#    print "kinfai"
    if line[0:5]=="track":
        print line.strip() + "\tuseScore=1"
        continue
    ls = line.strip().split('\t')
    gene_name = ls[1].split('.')[0]
    if not score_dt.has_key(gene_name):
         print "find no expression for GENE", gene_name
         output.write( line )
         continue

    if not score_dt[gene_name].has_key(ls[3]):
         print "find no expression for TRANSCRIPT",ls[1]
         output.write( line )
         continue

    ls[4] = str(score_dt[gene_name][ls[3]])
    ls[3] = ls[3]+"|"+ls[4]
    ls[8] = col_str
    for index in index_ls:
        ls[8] = design_col(int(ls[4]),0,255,ls[8],index)

    output.write( '\t'.join(ls) + '\n' )
gpd_file.close()
output.close()
