#!/usr/bin/python

import sys
import os
from numpy import *

if len(sys.argv) >= 4:
    refFlat_filename = sys.argv[1]
    junction_filename = sys.argv[2]
    exon_construction_junction_span = int(sys.argv[3])
else:
    print("usage: python exon_construction.py refFlat.txt junction.bed")
    print("or ./exon_construction.py refFlat.txt junction.bed")
    sys.exit(1)
################################################################################
def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    return path, filename

def check_jun(jun_list1,jun_list2):
    if (not jun_list1[1]<jun_list2[0]) and (not jun_list2[1]<jun_list1[0]):
        return 1
    return 0

################################################################################
def merge_generegion(d):
    result = {}
    ls = d.keys()
    ls.sort()
    ref_region = [ ls[0],d[ls[0]] ]
    for start in ls:
        end=int(d[start])
        start = int(start)

        if check_jun(ref_region,[start,end]) == 0:
            result[ref_region[0]]=ref_region[1]
            ref_region[0] = start
            ref_region[1] = end
        else:
            ref_region[0] = min(ref_region[0],start)
            ref_region[1] = max(ref_region[1],end)
    result[ref_region[0]]=ref_region[1]
    return result
    
################################################################################
def extend_generegion(d, chr_name, jun_dt):
    for jun_start in jun_dt[chr_name]:

        if not d.has_key(jun_start):
            d[jun_start]=max((jun_dt[chr_name][jun_start]))
        else:
            d[jun_start] = max(d[jun_start], max((jun_dt[chr_name][jun_start])) )
    return merge_generegion(d)
################################################################################

def print_jun(chr_name,jun_start, jun_end):
    output_ls = [chr_name,str(jun_start-50),str(jun_end+50),"novel_long_jun","50","+",str(jun_start-50),str(jun_end+50),"0,0,0","2","50,50"]
    jun_L = jun_end - jun_start +50
    output_ls.append("0,"+str(jun_L))
    print  "\t".join(output_ls)

################################################################################

ref_gene_region_dt={}
ref_jun_dt = {}
ref=open(refFlat_filename,'r')

for refline in ref:
    refline_list=refline.strip().split("\t")
    exon_start_list=refline_list[9].strip(',').split(',')
    exon_end_list=refline_list[10].strip(',').split(',')
    start = int(refline_list[4])
    end = int(refline_list[5])
    chr_name = refline_list[2]
    ID = refline_list[1]
    gene_name = refline_list[0]
    strand=refline_list[3]

    i = 1 
    for jun_start in exon_end_list[:-1]:
        locus = chr_name + ":" + jun_start + "-" + exon_start_list[i]
        ref_jun_dt[locus]=gene_name
        i += 1
    
    if not ref_gene_region_dt.has_key(chr_name):
        ref_gene_region_dt[chr_name] = {}
    if not ref_gene_region_dt[chr_name].has_key(gene_name):
        ref_gene_region_dt[chr_name][gene_name] = {}
    if not ref_gene_region_dt[chr_name][gene_name].has_key(start):
        ref_gene_region_dt[chr_name][gene_name][start] = end
    else:
        if ref_gene_region_dt[chr_name][gene_name][start] < end:
            ref_gene_region_dt[chr_name][gene_name][start] = end

ref.close()
####Define ref gene boundary########################################################

ref_gene_exon_range = {}
ref_gene_exon_range_locus_key = {}
for chr_name in ref_gene_region_dt:
    ref_gene_exon_range[chr_name] = {}
    ref_gene_exon_range_locus_key[chr_name] = {}
    for gene_name in ref_gene_region_dt[chr_name]:
        ref_gene_exon_range[chr_name][gene_name] = merge_generegion(ref_gene_region_dt[chr_name][gene_name])

        gene_start = ref_gene_exon_range[chr_name][gene_name].keys()[0]
        gene_end = ref_gene_exon_range[chr_name][gene_name][gene_start]
        ref_gene_exon_range_locus_key[chr_name][gene_start] = {}
        ref_gene_exon_range_locus_key[chr_name][gene_start][gene_end] = gene_name

################################################################################

jun_dt = {}
novel_jun_dt = {}
jun_file = open(junction_filename,'r')
for line in jun_file:
    if line[0:5]=="track":
        continue
    ls = line.strip().split('\t')

    chr_name = ls[0]
    thickness=ls[10].split(',')
    leftpos=int(ls[1])+int(thickness[0])
    rightpos=int(ls[2])-int(thickness[1])

    name = ls[3]  
    num_indep=int(name.split('_')[1].split(']')[0])
    num_uniq=int(name.split('](')[1].split('/')[0])

#    print leftpos, rightpos, num_indep, num_uniq

    if not (num_indep > 1 and num_uniq > 0):
        continue

    locus = chr_name + ":" + str(leftpos) + "-" + str(rightpos)
    if not ref_jun_dt.has_key(locus):
        if not novel_jun_dt.has_key(chr_name):
            novel_jun_dt[chr_name] = {}
        if not novel_jun_dt[chr_name].has_key(leftpos):
            novel_jun_dt[chr_name][leftpos] = []
        novel_jun_dt[chr_name][leftpos].append(rightpos)
    else:
        if not jun_dt.has_key(chr_name):
            jun_dt[chr_name] = {}
        if not jun_dt[chr_name].has_key(leftpos):
            jun_dt[chr_name][leftpos] = []
        jun_dt[chr_name][leftpos].append(rightpos)
jun_file.close()

#######add filtered novel jun##################################################

for chr_name in novel_jun_dt:
    jun_start_ls = novel_jun_dt[chr_name].keys()
    if not ref_gene_exon_range_locus_key.has_key(chr_name):
        continue
    ref_gene_exon_range_start_ls = ref_gene_exon_range_locus_key[chr_name].keys()
    ref_gene_exon_range_start_ls.sort()
    L_ref_gene_exon_range_start_ls = len(ref_gene_exon_range_start_ls)
    jun_start_index_ls = searchsorted( ref_gene_exon_range_start_ls, jun_start_ls )
    i = 0 
    for jun_start_index in jun_start_index_ls:
        jun_start = jun_start_ls[i]
        jun_end_ls = novel_jun_dt[chr_name][jun_start]
        for jun_end in jun_end_ls:


            j = 0
            while jun_start_index+j < L_ref_gene_exon_range_start_ls and  jun_end > ref_gene_exon_range_start_ls[jun_start_index+j]:
                j+=1

            if check_jun([jun_start,jun_end],[ ref_gene_exon_range_start_ls[jun_start_index-1], ref_gene_exon_range_locus_key[chr_name][ref_gene_exon_range_start_ls[jun_start_index-1]] ]) == 1:
                j+=1

            if j <= exon_construction_junction_span:
                if not jun_dt.has_key(chr_name):
                    jun_dt[chr_name] = {}
                if not jun_dt[chr_name].has_key(jun_start):
                    jun_dt[chr_name][jun_start] = []
                jun_dt[chr_name][jun_start].append(jun_end)
            else:
                print_jun(chr_name,jun_start, jun_end)
        i+=1
# sys.exit(1)

################################################################################
ref=open(refFlat_filename,'r')
exon_start_end_dt = {}
jun_start_end_dt = {}

for refline in ref:
    refline_list=refline.strip().split("\t")
    exon_start_list=refline_list[9].strip(',').split(',')
    exon_end_list=refline_list[10].strip(',').split(',')
    chr_name = refline_list[2]
    gene_name = refline_list[1]

    if not exon_start_end_dt.has_key(chr_name):
        exon_start_end_dt[chr_name]={}
    if not exon_start_end_dt[chr_name].has_key(gene_name):
        exon_start_end_dt[chr_name][gene_name]=[[],[]]

    exon_start_end_dt[chr_name][gene_name][0].append(int(exon_start_list[0]))
    exon_start_end_dt[chr_name][gene_name][1].append(int(exon_end_list[-1]))

    ######################################

    if not jun_start_end_dt.has_key(chr_name):
        jun_start_end_dt[chr_name]={}
    if not jun_start_end_dt[chr_name].has_key(gene_name):
        jun_start_end_dt[chr_name][gene_name]=[[],[]]

    if len(exon_end_list)>1:
        jun_start_end_dt[chr_name][gene_name][0].append(int(exon_end_list[0]))
        jun_start_end_dt[chr_name][gene_name][1].append(int(exon_start_list[-1]))
    else:
        jun_start_end_dt[chr_name][gene_name][0].append(int(exon_start_list[0]))
        jun_start_end_dt[chr_name][gene_name][1].append(int(exon_end_list[-1]))
          
ref.close()

####Define gene boundary########################################################

gene_exon_range = {}
gene_region = {}
for chr_name in exon_start_end_dt:
    gene_exon_range[chr_name] = {}
    for gene_name in exon_start_end_dt[chr_name]:

        k = 0 
        for start in exon_start_end_dt[chr_name][gene_name][0]:
            end = exon_start_end_dt[chr_name][gene_name][1][k]
            k += 1
            if not gene_exon_range[chr_name].has_key(start):
                gene_exon_range[chr_name][start] = end
            else:
                if gene_exon_range[chr_name][start] < end:
                    gene_exon_range[chr_name][start] = end

    gene_region[chr_name] = merge_generegion(gene_exon_range[chr_name])
    if jun_dt.has_key(chr_name):
        gene_region[chr_name] = extend_generegion(gene_region[chr_name], chr_name, jun_dt)

####Add junctions from RefSeq###################################################

ref_gene_boundary_dt = {}

ref=open(refFlat_filename,'r')
for refline in ref:
    refline_list=refline.split()
    Exon_start=int(refline_list[4])
    Exon_end=int(refline_list[5])
    Exon_start_list=refline_list[9].strip(",").split(',')
    Exon_end_list=refline_list[10].strip(",").split(',')
    strand=refline_list[3]
    chr_name = refline_list[2]
    i = 1 
    output_str = ""

    if not jun_dt.has_key(chr_name):
        jun_dt[chr_name] = {}

    for jun_start in Exon_end_list[:-1]:
        jun_start = int(jun_start)
        jun_end = int(Exon_start_list[i])
        if not jun_dt[chr_name].has_key(jun_start):
            jun_dt[chr_name][jun_start] = []
        jun_dt[chr_name][jun_start].append(jun_end)
        i += 1

    if not jun_dt[chr_name].has_key (Exon_start):
        jun_dt[chr_name][Exon_start] = []
    jun_dt[chr_name][Exon_start].append(Exon_end)

################################################################################
exon_dt = {}

for chr_name in jun_dt:
    if not gene_region.has_key(chr_name):
        print "not gene region", chr_name
        continue
    region_start_ls = gene_region[chr_name].keys()
    region_start_ls.sort()
    i = 0
    gene_region_start = region_start_ls[i]
    gene_region_end = gene_region[chr_name][gene_region_start]

    jun_start_ls = jun_dt[chr_name].keys()
    jun_start_ls.sort()

    exon_dt[chr_name]={}
    temp_exon_boundary = []

    for jun_start in jun_start_ls:

        while jun_start >= gene_region_end:
            if i < len(region_start_ls):
                locus = chr_name + ':' + str(gene_region_start) + '-' + str(gene_region_end)
                exon_dt[chr_name][locus]=temp_exon_boundary 
                temp_exon_boundary = []
                i += 1
                if i < len(region_start_ls):
                    gene_region_start = region_start_ls[i]
                    gene_region_end = gene_region[chr_name][gene_region_start]

        jun_end_ls = jun_dt[chr_name][jun_start]
        jun_end_ls.sort()



        if jun_end_ls[-1] < gene_region_start:
            print "skip", chr_name, ":" , region_start_ls[i] , "-" , gene_region[chr_name][region_start_ls[i]]
            continue

        for jun_end in jun_end_ls:
            if jun_start >= gene_region_start and jun_end <= gene_region_end:
                temp_exon_boundary.append(jun_start)
                temp_exon_boundary.append(jun_end)

    if i < len(region_start_ls):
        locus = chr_name + ':' + str(gene_region_start) + '-' + str(gene_region_end)
        exon_dt[chr_name][locus]=temp_exon_boundary 

#        print locus, exon_dt[chr_name][locus]

#    while  i == len(region_start_ls)-1: 
#        print chr_name , ":" , region_start_ls[i] , "-" , gene_region[chr_name][region_start_ls[i]]
#        i+=1
refFlat_path,refFlat_filename=GetPathAndName(refFlat_filename)

exon_output = open( junction_filename + "_" + refFlat_filename+".exon",'w')
for chr_name in exon_dt:
    for locus in exon_dt[chr_name]:

        output_cmd = locus

        locus_ls = locus.split(':')
        locus_ls_ls = locus_ls[1].split('-')
        start = int(locus_ls_ls[0])
        end = int(locus_ls_ls[1])
        chr_name = locus_ls[0]


        pt_ls = list(set(exon_dt[chr_name][locus]))
        pt_ls = list(set(pt_ls))
        pt_ls.sort()
#        print len(pt_ls)
        for pt in pt_ls:
            output_cmd = output_cmd + '\t' + str(pt) 
        exon_output.write(output_cmd+ '\n')
exon_output.close()




