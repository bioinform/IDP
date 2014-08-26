#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 4:
    pos_filename = sys.argv[1]
    neg_filename = sys.argv[2]
    FDR = float(sys.argv[3])
    input_filename = sys.argv[4]
    output_filename = sys.argv[5]
else:
    print("usage: ./ROC.py pos_filename neg_filename FDR output_filename")
    print("or ")
    sys.exit(1)
################################################################################
#chr9:127997126-128003666.1      583.3308        583.3308

################################################################################
neg_perc_ls = []
neg_file = open(neg_filename,'r')
for line in neg_file:
    ls = line.strip().split("\t")
    ID = ls[0]
    iso_exp = float(ls[1])
    gene_exp = float(ls[2])
    if (gene_exp > 0):
        neg_perc_ls.append(iso_exp/gene_exp)
    else:
        neg_perc_ls.append(0)
neg_file.close()

neg_perc_ls.sort()
neg_perc_ls.reverse()
n = int(FDR*len(neg_perc_ls))
threshold = neg_perc_ls[n]
print "negative:",n,len(neg_perc_ls),float(n)/len(neg_perc_ls)
print "threshold:",threshold
################################################################################
output = open(output_filename,'w')
input_file = open(input_filename,'r')
i=0
I=0
for line in input_file:
    I+=1
    ls = line.strip().split("\t")
    ID = ls[0]
    iso_exp = float(ls[1])
    gene_exp = float(ls[2])
    if gene_exp > 0 and iso_exp/gene_exp > threshold:
        output.write(line)
        i+=1
input_file.close()
print "input:",i,I,float(i)/I
#print "FDR:", float(n)/(n+i)
output.close()
################################################################################
