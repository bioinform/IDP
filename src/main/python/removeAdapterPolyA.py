#!/usr/bin/python

import sys
import os

if len(sys.argv)>=6:
    in_filename = sys.argv[1]
    idx_filename =sys.argv[2]
    Three_filename =sys.argv[3]
    Five_filename = sys.argv[4]
    polyA_min_len = int(sys.argv[5])
    out_filename = sys.argv[6]
     
else:
    print("usage: removeAdapterPolyA.py in_filename idx_filename Three_filename Five_filename polyA_min_len out_filename")
    sys.exit(1)

################################################################################
primer_margin = 11

################################################################################

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
Five_dt={}
Five = open(Five_filename,'r')
# Read the header
Five.readline()
for line in Five:
    ls=line.strip().split('\t')
    readname = ls[0]
    if not len_dict.has_key(readname):
        continue
    pos = int(ls[1]) - 1 #seqmap is 1-based
    strand = ls[-1]
    L = len(ls[2].replace('_','') )
    ls[2] = L
    primer_L = len(ls[4])

    if not Five_dt.has_key(readname):
        Five_dt[readname]=[]

    if ((strand == '+') and (pos <= primer_margin)):
        Five_dt[readname].append(ls)
    if ((strand == '-') and (pos >= (len_dict[readname]- primer_margin - primer_L))):
        Five_dt[readname].append(ls)

Five.close()

################################################################################

Three_dt={}
Three= open(Three_filename,'r')
# Read the header
Three.readline()

for line in Three:
    ls=line.strip().split('\t')
    readname = ls[0]
    if not len_dict.has_key(readname):
        continue
    pos = int(ls[1]) - 1  #seqmap is 1-based
    strand = ls[-1]
    L = len( ls[2].replace('_','') )
    ls[2] = L
    primer_L = len(ls[4])

    if not Three_dt.has_key(readname):
        Three_dt[readname]=[]

    if ((strand == '-') and (pos <= primer_margin)) :
        polyT_pos = pos + L
        NT = 0
        if idx_dt[readname].has_key(polyT_pos):
            NT = idx_dt[readname][polyT_pos]
        if ((cps_dt[readname][polyT_pos] == 'T') and 
            (NT >= polyA_min_len)):
            Three_dt[readname].append(ls)

    if ((strand == '+') and (pos >= (len_dict[readname]- primer_margin - primer_L))):
        polyA_pos = pos - 1
        NA = 0
        if idx_dt[readname].has_key(polyA_pos):
            NA = idx_dt[readname][pos-1]
        if ((cps_dt[readname][polyA_pos] =='A') and (NA >= polyA_min_len)):
            Three_dt[readname].append(ls)

Three.close()

################################################################################

#Output: strand and 1-based position after/before primer sequence
def cutmost_five(temp_ls_ls, FL):
    s = ''
    ref_cut_L = 0
    for temp_ls in temp_ls_ls:
        pos = int(temp_ls[1])
        if temp_ls[-1]== '+':
            cut_L = pos + temp_ls[2] - 1
            temp_keep_pt = pos + temp_ls[2]
        else:
            cut_L = FL - pos  + 1
            temp_keep_pt = pos - 1 
        if ref_cut_L <= cut_L:
            keep_pt = temp_keep_pt
            s = temp_ls[-1]
    return s, keep_pt
 
#Output: strand and selected primer mapping            
def cutmost_three(temp_ls_ls, FL):
    s = ''
    ref_cut_L = 0
    result_ls = []
    for temp_ls in temp_ls_ls:
        pos = int(temp_ls[1])
        if temp_ls[-1] == '+':
            cut_L = FL - pos + 1
        else:
            cut_L = pos + temp_ls[2] - 1
  
        
        if ref_cut_L <= cut_L:
            s = temp_ls[-1]
            result_ls = temp_ls 
    return s, result_ls

def process(Three_ls_ls, Five_ls_ls, FL, readname):
    if len(Three_ls_ls) == 0:
        strand = ''
    
    elif len(Three_ls_ls) == 1:
        strand = Three_ls_ls[0][-1]

    else:
        strand, result_ls = cutmost_three(Three_ls_ls, FL)
        Three_ls_ls = [result_ls]  # Update it with the selected entry
        print ">1 polyA/T", readname
#        return 1,FL

    # No PolyA is detected
    if strand == '':
        if len(Five_ls_ls) == 0:
            return 1,FL,strand
        
        strand, keep_pt = cutmost_five(Five_ls_ls,len_dict[readname])        
        if strand =='-':
            return 1, keep_pt, ''
        else:
            return keep_pt, FL, ''

    else:
        temp_five_ls_ls = []
        for temp_ls in Five_ls_ls:
            if temp_ls[-1] == strand:
                temp_five_ls_ls.append(temp_ls)

        if strand == '+':
            right_keep_pt = int(Three_ls_ls[0][1]) - 2
            if len(temp_five_ls_ls) == 0:
                return 1, right_keep_pt, strand
            else:
                strand, keep_pt = cutmost_five(temp_five_ls_ls, FL)
                return keep_pt,right_keep_pt,strand
        else:
            left_keep_pt = int(Three_ls_ls[0][1]) + Three_ls_ls[0][2] + 1
            if len(temp_five_ls_ls)==0:
                return left_keep_pt, FL, strand
            else:
                strand, keep_pt = cutmost_five(temp_five_ls_ls, FL)
                return left_keep_pt,keep_pt,strand

################################################################################
def makelist(seq, temp_idx_dt):
    result_list = list(seq)   
    for item in temp_idx_dt:
        result_list[item]=result_list[item]*temp_idx_dt[item]
    return result_list
################################################################################

infile = open(in_filename,'r')
idx_file =open(idx_filename,'r')

outfile = open(out_filename,'w')
threeend_file = open(out_filename+'.3','w')
for line in infile:
    
    if line[0]=='>':
         outfile.write(line)
         readname=line.strip().strip('>')
    else:
         if not Three_dt.has_key(readname):
             Three_dt[readname]=[]
         if not Five_dt.has_key(readname):
             Five_dt[readname]=[]
         start, end, strand = process(Three_dt[readname], Five_dt[readname], len_dict[readname], readname)

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
    
