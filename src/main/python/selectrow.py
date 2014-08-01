#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 4 :
    filename = sys.argv[1]
    row_i = int(sys.argv[2])-1
    target_ls_filename = sys.argv[3]
    output_filename = sys.argv[4]
else:
    print("usage: python selectrow.py filename row_i target_ls_filename")
    print("or ./selectrow.py filename row_i target_ls_filename")
    sys.exit(1)
################################################################################
file = open(filename,'r')
dt = {}
for line in file:
    ls=line.strip().split('\t')
    if not dt.has_key(ls[row_i]):
        dt[ ls[row_i] ] = []
    dt[ ls[row_i] ].append( line.strip() )
file.close()

################################################################################
output = open(output_filename,'w')
target_ls_file = open(target_ls_filename, 'r')
for line in target_ls_file:
    id = line.strip()
    if not dt.has_key(id):
        print id
        continue
    if len(dt[id])>1:
        print id + '\t' + str(len(dt[id]))
    for item in dt[id]:
        output.write( item + '\n')
output.close()
target_ls_file.close()
