#!/usr/bin/python

import sys
import os

if len(sys.argv)>=6:
    in_filename = sys.argv[1]
    idx_filename =sys.argv[2]
    ACTGTACAGT_filename =sys.argv[3]
    ACACTCTGT_filename =sys.argv[4]
    out_filename = sys.argv[5]

else:
    print("usage:./removeAdapterPolyA.py in_filename idx_filename ACTGTACAGT_filename ACACTCTGT_filename out_filename")
    print("or ./removeAdapterPolyA.py test_LR.fa.cps test_LR.fa.idx test_LR.fa.cps_ACTGTACAGT.out test_LR.fa.cps_ACACTCTGT.out LR_nopolyA.fa")
    sys.exit(1)
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
ACTGTACAGT_dt={}
ACTGTACAGT = open(ACTGTACAGT_filename,'r')
i=0
for line in ACTGTACAGT:
    if i==0:
        i+=1
        continue
    ls=line.strip().split('\t')
    readname = ls[0]
    pos = int(ls[1])
    strand = ls[-1]
    L =len( ls[2].replace('_','') )
    ls[2]=L

    if not ACTGTACAGT_dt.has_key(readname):
        ACTGTACAGT_dt[readname]=[]

    if strand == '+' and pos<=11 :
        ACTGTACAGT_dt[readname].append(ls)
    if strand == '-' and pos>=len_dict[readname]-19:
        ACTGTACAGT_dt[readname].append(ls)

ACTGTACAGT.close()

################################################################################

ACACTCTGT_dt={}
ACACTCTGT= open(ACACTCTGT_filename,'r')
i=0
for line in ACACTCTGT:
    if i==0:
        i+=1
        continue
    ls=line.strip().split('\t')
    readname = ls[0]
    pos = int(ls[1])
    ls[1]=pos
    strand = ls[-1]
    L =len( ls[2].replace('_','') )
    ls[2]=L

    if not ACACTCTGT_dt.has_key(readname):
        ACACTCTGT_dt[readname]=[]

    if strand == '-' and pos<=16 :
        polyT_pos = pos+L-2
        NT = 1
        if idx_dt[readname].has_key(polyT_pos):
            NT = idx_dt[readname][polyT_pos]
        if cps_dt[readname][polyT_pos] =='T' and NT>=10:
            ACACTCTGT_dt[readname].append(ls)

    if strand == '+' and pos>=len_dict[readname]-23:
        NA = 1
        if idx_dt[readname].has_key(pos-1):
            NA = idx_dt[readname][pos-1]
        if cps_dt[readname][pos-1] =='A' and NA>=10:
            ACACTCTGT_dt[readname].append(ls)

ACACTCTGT.close()

################################################################################

def cutmost(temp_ls_ls,FL):
    s=''
    ref_cut_L = 0
    for temp_ls in temp_ls_ls:
        pos = int(temp_ls[1])
        if temp_ls[-1]=='+':
            cut_L = pos-1 + temp_ls[2]
            temp_keep_pt = pos + temp_ls[2]
        else:
            cut_L = FL - pos +1 
            temp_keep_pt = pos-1
        if ref_cut_L<cut_L:
            keep_pt = temp_keep_pt
            s = temp_ls[-1]
    return s, keep_pt
            
def cutmost2(temp_ls_ls,FL):
    s=''
    ref_cut_L = 0
    result_ls = []
    for temp_ls in temp_ls_ls:
        pos = int(temp_ls[1])
        if temp_ls[-1]=='+':
            cut_L = FL - pos +1 
            temp_keep_pt = pos-1
        else:
            cut_L = pos-1 + temp_ls[2]
            temp_keep_pt = pos + temp_ls[2]
        if ref_cut_L<cut_L:
            keep_pt = temp_keep_pt
            s = temp_ls[-1]
            result_ls = temp_ls 
    return s, result_ls

def process(ACACTCTGT_ls_ls,ACTGTACAGT_ls_ls,FL,readname):
    if len(ACACTCTGT_ls_ls)==0:
        strand = ''
    elif len(ACACTCTGT_ls_ls)==1:
        strand = ACACTCTGT_ls_ls[0][-1]

    else:
        strand, result_ls = cutmost2(ACACTCTGT_ls_ls,FL)
        ACACTCTGT_ls_ls = [result_ls]
        print ">1 polyA/T", readname
#        return 1,FL

    if strand == '':
        if len(ACTGTACAGT_ls_ls)==0:
            return 1,FL,strand
        strand, keep_pt = cutmost(ACTGTACAGT_ls_ls,len_dict[readname])        
        if strand =='-':
            return 1,keep_pt,''
        else:
            return keep_pt,FL,''
    else:
        temp_ls_ls = []
        for temp_ls in ACTGTACAGT_ls_ls:
            if temp_ls[-1]==strand:
                temp_ls_ls.append(temp_ls)
        if strand == '+':
            right_keep_pt = ACACTCTGT_ls_ls[0][1]-1
            if len(temp_ls_ls)==0:
                return 1,right_keep_pt,strand
            else:
                strand, keep_pt = cutmost(temp_ls_ls,FL)
                return keep_pt,right_keep_pt,strand
        else:
            left_keep_pt = ACACTCTGT_ls_ls[0][1]+ACACTCTGT_ls_ls[0][2]
            if len(temp_ls_ls)==0:
                return left_keep_pt,FL,strand
            else:
                strand, keep_pt = cutmost(temp_ls_ls,FL)
                return left_keep_pt,keep_pt,strand


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
         if not ACACTCTGT_dt.has_key(readname):
             ACACTCTGT_dt[readname]=[]
         if not ACTGTACAGT_dt.has_key(readname):
             ACTGTACAGT_dt[readname]=[]
         start,end,strand = process(ACACTCTGT_dt[readname],ACTGTACAGT_dt[readname],len_dict[readname],readname)

# strand == "+", end is polyAT
# strand == "-", start is polyAT
# strand == "",  no polyAT
         if strand != "":
             threeend_file.write(readname + "\t" + strand + "\n" )

         cps_list = makelist(line.strip(), idx_dt[readname])
         final_seq = ''.join(cps_list[start-1:end])
         outfile.write(final_seq+'\n')

outfile.close()
infile.close()
idx_file.close()
threeend_file.close()

################################################################################
    
