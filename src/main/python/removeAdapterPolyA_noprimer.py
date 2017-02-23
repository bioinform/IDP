#!/usr/bin/python

import sys
import os

if len(sys.argv)>=4:
    in_filename = sys.argv[1]
    idx_filename =sys.argv[2]
    polyA_min_len = int(sys.argv[3])
    out_filename = sys.argv[4]
     
else:
    print("usage: removeAdapterPolyA.py in_filename idx_filename Three_filename Five_filename polyA_min_len out_filename")
    sys.exit(1)

################################################################################
primer_margin = 11

################################################################################
print("Parse idx file")
idx_dt={}
idx_file =open(idx_filename,'r')
for line in idx_file:
    ls=line.strip().split('\t')
    readname =ls[0]
    idx_dt[readname]={}
    if len(ls)>1:
        pos_ls=ls[1].split(',')
        n_ls = ls[2].split(',')
        i=0
        for pos in pos_ls:
            idx_dt[readname][int(pos)]=int(n_ls[i])
            i+=1
idx_file.close()

################################################################################
print("Parse cps file")
cps_dt={}
len_dict={}
infile = open(in_filename,'r')
for line in infile:
    if line[0] == ">":
        readname = line.strip().strip('>')
        if not len_dict.has_key(readname):
            len_dict[readname]=0
            cps_dt[readname]=''
        else:
            print "Error: repeat readnames:", readname
        
    else:
        len_dict[readname] += len(line.strip())
        cps_dt[readname]=cps_dt[readname]+line.strip()
infile.close()


################################################################################
def process(FL, readname):

    # Check polyA or polyT seq at either end of the read
    if idx_dt[readname].has_key(FL - 1):
        NA = idx_dt[readname][FL - 1]
        if ((cps_dt[readname][FL - 1] == 'A') and 
            (NA >= polyA_min_len)):
            return 1, FL - 1, '+'
    
    if idx_dt[readname].has_key(0):
        NT = idx_dt[readname][0]
        if ((cps_dt[readname][0] == 'T') and 
            (NT >= polyA_min_len)):
            return 2, FL, '-'
    
    return 1, FL, ''

################################################################################
def makelist(seq, temp_idx_dt):
    result_list = list(seq)   
    for item in temp_idx_dt:
        result_list[item]=result_list[item]*temp_idx_dt[item]
    return result_list
################################################################################
print("Generate output file")
infile = open(in_filename,'r')
idx_file =open(idx_filename,'r')

outfile = open(out_filename,'w')
threeend_file = open(out_filename+'.3','w')
for line in infile:
    
    if line[0]=='>':
         outfile.write(line)
         readname=line.strip().strip('>')
    else:
         
         start, end, strand = process(len_dict[readname], readname)

# strand == "+", end is polyAT
# strand == "-", start is polyAT
# strand == "",  no polyAT
         if strand != "":
             threeend_file.write(readname + "\t" + strand + "\n" )

         cps_list = makelist(line.strip(), idx_dt[readname])
         final_seq = ''.join(cps_list[start - 1: end])
         outfile.write(final_seq + '\n')

outfile.close()
infile.close()
idx_file.close()
threeend_file.close()

################################################################################
    
