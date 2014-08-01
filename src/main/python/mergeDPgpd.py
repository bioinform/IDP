#!/usr/bin/python

import sys
import os
from numpy import *
################################################################################

if len(sys.argv) >= 4:
    name_filename = sys.argv[1]
    gpd_filename =  sys.argv[2]
    output_filename =  sys.argv[3]
else:
    print("usage: python name_filename gpd_filename output_filename")
    print("or ./g name_filename gpd_filename  output_filename")
    sys.exit(1)

################################################################################

output = open(output_filename,"w")
name_dt={}
name_file = open(name_filename,"r")
for line in name_file:
    ls = line.strip().split("\t")
    gene = ls[1]
    ID = ls[0]
    refgene = ls[2]
    refID = ls[3]
    if refID == "-":
        name_dt[ID] = ls[1:]
name_file.close()

tag = open(gpd_filename,'r')
for line in tag:
    ls = line.strip().split('\t')
    gene = ls[0]
    ID = ls[1]
    chr_name = ls[2]
    if ls[8]=="1":
        continue
    if not name_dt.has_key(ID):
        continue
    if name_dt[ID][1]!="-":
        ls[0] = name_dt[ID][1]
    output.write( '\t'.join(ls) + "\n" )
tag.close()
output.close()

################################################################################
















