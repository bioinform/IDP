#!/usr/bin/python

import sys
import os
from numpy import *
from scipy import stats

if len(sys.argv) >= 4 :
    ref_filename = sys.argv[1]
    tag_filename =sys.argv[2]
    Npt = int(sys.argv[3])
    Nbin = int(sys.argv[4])
else:
    print("usage: ~/3seq/bin/exp_len_density.py multiexon_refFlat.txt_positive_known_intact_SM.fa.bestpsl.gpd_refFlat.txt_exp_len multiexon_refFlat.txt_positive_known_intact_SM.fa.bestpsl.gpd_refFlat.txt_exp_len 100")
    print("or ")
    sys.exit(1)
################################################################################

ref = open(ref_filename,'r')
len_dt = {}
for line in ref:
    ls = line.strip().split("\t")
    L = int(ls[2])
    if not len_dt.has_key(L):
        len_dt[L]=[]
    len_dt[L].append(ls)
ref.close()

################################################################################

def getdensity(len_ls,len_dt, L,Npt):
    result= [] 
    index = searchsorted(len_ls,L,side='right')

    left_index = index - 1 
    right_index = index 
    left_L = len_ls[left_index]
    right_L = len_ls[right_index]

    r_left_L = L - left_L
    r_right_L = right_L - L

    left_iso_ls = []
    right_iso_ls = []

    if left_L > smallnum:
        left_iso_ls =len_dt[left_L]
        
    if right_L < largenum:
        right_iso_ls = len_dt[right_L]

    len_left_iso_ls = len(left_iso_ls)
    len_right_iso_ls = len(right_iso_ls)

    if len_left_iso_ls + len_right_iso_ls > Npt:
        if r_left_L < r_right_L:
            if len_left_iso_ls > Npt:
                return left_iso_ls[:Npt]
            else:
                result.extend(left_iso_ls)
                result.extend(right_iso_ls[:Npt-len_left_iso_ls])
                return result
        else:
            if len_right_iso_ls > Npt:
                return right_iso_ls[:Npt]
            else:
                result.extend(right_iso_ls)
                result.extend(left_iso_ls[:Npt-len_right_iso_ls])
                return result

    n = len(result)

    while len(result)<Npt:

        if r_left_L < r_right_L:
            while r_left_L < r_right_L and len(result)<Npt:
                result.extend(left_iso_ls)
                left_index -= 1
                left_L = len_ls[left_index]
                if left_L > smallnum:
                    left_iso_ls =len_dt[left_L]
                r_left_L = L - left_L

        else:
            while r_left_L >= r_right_L and len(result)<Npt:
                result.extend(right_iso_ls)
                right_index += 1
                right_L = len_ls[right_index]
                if right_L < largenum:
                     right_iso_ls = len_dt[right_L]
                r_right_L = right_L - L
    return result[:Npt]

################################################################################
def calculate_b(pt,npt):
    RPKM_ls = []
    I_ls =[]
    L_ls = []
    for item in pt:
       RPKM_ls.append( float(item[3]) )
       I_ls.append(int(item[4]))
       L_ls.append(int(item[2]))

    temp_a = array([RPKM_ls,I_ls])

    temp_a =transpose(temp_a)
    temp_a_sorted = transpose( sorted(temp_a, key=lambda a_entry: a_entry[0]) )
    RPKM_med_ls = []
    D_rate_ls = []
    i = 0
    L_pt = len(pt)
    while i < L_pt:
        RPKM_med_ls.append( median( temp_a_sorted[0][i:i+npt] ) )
        D_rate_ls.append( 1-mean( temp_a_sorted[1][i:i+npt] ) )
        i += npt
    gradient, intercept, r_value, p_value, std_err = stats.linregress(RPKM_med_ls, D_rate_ls)
    return gradient, intercept, r_value, p_value, std_err,std(L_ls)

def printout(pt):
    result = []
    s = 0
    for item in pt:
        result.append(str(item[2]))
        s += float(item[2])
    print '\t'.join(result)
    ave2 = s/len(result)


    result = []
    s = 0
    for item in pt:
        result.append(str(item[3]))
        s += float(item[3])
    print '\t'.join(result)
    ave3 = s/len(result)

    return ave2, ave3

len_ls = len_dt.keys()
largenum = 1e10
smallnum = -1e10

len_ls.append(largenum)
len_ls.append(smallnum)

len_ls.sort()

L=0
while L<0:
    pt = getdensity(len_ls,len_dt, L,Npt)
    if len(pt)!=Npt:
        sys.exit(1)
    gradient, intercept, r_value, p_value, std_err,std_L = calculate_b(pt,Npt/Nbin)
    print '\t'.join([str(L),str(gradient), str(intercept), str(r_value), str(p_value), str(std_err),str(std_L)])
    L+=1

#sys.exit(1)

tag = open(tag_filename,'r')
for line in tag:
    ls = line.strip().split("\t")

    exon_start_list=ls[9].strip(',').split(',')
    exon_end_list=ls[10].strip(',').split(',')

    L = 0
    i=0
    for start in exon_start_list:
        start =int(start)
        end = int(exon_end_list[i])
        L += (end - start)
        i += 1

    pt = getdensity(len_ls,len_dt, L,Npt)
    if len(pt)!=Npt:
        sys.exit(1)
    gradient, intercept, r_value, p_value, std_err,std_L = calculate_b(pt,Npt/Nbin)
    print '\t'.join([str(L),str(gradient), str(intercept), str(r_value), str(p_value), str(std_err),str(std_L)])
tag.close()

