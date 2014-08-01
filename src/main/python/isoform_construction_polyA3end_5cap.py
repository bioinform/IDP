#!/usr/bin/python

from numpy import *
import sys
import os
from copy import *
exon_len_limit = 4000

if len(sys.argv) >= 11:
    refFlat_filename = sys.argv[1]
    junction_filename = sys.argv[2]
    exon_boundary_filename = sys.argv[3]
    CAGE_filename = sys.argv[4]
    nonredundant_LR_gpd_filename = sys.argv[5]
    I_ref5end_isoformconstruction = sys.argv[6]
    I_ref3end_isoformconstruction = sys.argv[7]
    limit = int(sys.argv[8])
    exon_len_limit = int(sys.argv[9])
    output_filename = sys.argv[10]
#    I_fix3end = sys.argv[9]
else:
    print("usage: python isoform_construction.py refFlat.txt junction.bed exon_boundary_filename LR_gpd_filename")
    print("or /home/kinfai/3seq/bin/isoform_construction_polyA3end_5cap.py /home/3seq/kinfai/isoform_contrusction/refFlat.txt  /home/3seq/H1_SR_junction/junction_color.bed /home/kinfai/3seq/isoform_contrusction/junction_color.bed_refseq.exon LR_gpd_filename")
    sys.exit(1)

if CAGE_filename == "kinfai" and I_ref5end_isoformconstruction != "1":
    print "no CAGE data, you need reference 5 end"
    sys.exit(1)
#/home/kinfai/3seq/bin/isobin/isoform_construction_polyA3end_5cap.py /home/kinfai/temp/ref.gpd /home/kinfai/temp/SR.bed /home/kinfai/temp/SR.bed_ref.gpd.exon kinfai /home/kinfai/temp/junfil_compatible_LR_polyA3end.gpd 1 1 20 4000 > /home/kinfai/temp/isoform_construction.log

I_stop = 1

################################################################################

def check_3end(readname):
#       - is left
#       + is right
#       "" is no polyA3end
    I = readname.find('-')
    if I == -1:
        I = readname.find('+')
        if I == -1:
            return "",0
        else:
            return '+',int(readname[I+1:])
    else:
        return '-',int(readname[I+1:])
################################################################################

#gene_name	n
#a1-b1	a2-b2	a3-b3
#n1	n2	n3

def parse_5cap_file(CAGE_filename):
    if CAGE_filename == "kinfai":
        return {}
    result = {}
    CAGE_file = open(CAGE_filename,'r')
    i = 0
    for line in CAGE_file:
        ls = line.strip().split('\t')
        if ls[0].find(":") != -1 and ls[0].find("-") != -1:
            i = 0
        if 1:
            if i == 0:
                ls = line.strip().split('\t')
                gene_name = ls[0]
                chr_name = ls[1]
                N5cap = int(ls[2])
            elif i == 1:
                interval_str_ls = line.strip().split('\t')
            elif i == 2:
                Npeak_str_ls = line.strip().split('\t')
                j = 0
                for interval_str in interval_str_ls:
                    interval_ls = interval_str.split('-')
                    start = int(interval_ls[0])
                    end = int(interval_ls[1])
                    if not result.has_key(gene_name):
                         result[gene_name] = {}
                    result[gene_name][start]=[end,int(Npeak_str_ls[j])]
                    j += 1
            i +=1
    CAGE_file.close()
    return result

def define_5end(cap_dt,Nthreshold):
    fwd_5end_dt = {}
    rev_5end_dt = {}
    for gene_name in cap_dt:
        fwd_5end_dt[gene_name]=[]
        rev_5end_dt[gene_name]=[]
        for start in cap_dt[gene_name]:
            if cap_dt[gene_name][start][1] < Nthreshold:
                continue
            fwd_5end_dt[gene_name].append(start)
            rev_5end_dt[gene_name].append(cap_dt[gene_name][start][0])
    return fwd_5end_dt,rev_5end_dt
     
################################################################################

jun_start_dt = {}
jun_end_dt = {}
jun_strand_dt = {}

jun_file = open(junction_filename,'r')
for line in jun_file:
    if line[0:5]=="track":
        continue
    ls = line.strip().split('\t')

    chr_name = ls[0]
    strand = ls[5]
    thickness=ls[10].split(',')
    leftpos=int(ls[1])+int(thickness[0])
    rightpos=int(ls[2])-int(thickness[1])

    name = ls[3]  
    num_indep=int(name.split('_')[1].split(']')[0])
    num_uniq=int(name.split('](')[1].split('/')[0])

    locus = chr_name + ':' + str(leftpos) + '-' +str(rightpos)
    if not jun_strand_dt.has_key(locus):
        jun_strand_dt[locus] = strand

    if num_indep < 2 or num_uniq == 0:
        continue

    if not jun_start_dt.has_key(chr_name):
        jun_start_dt[chr_name] = {}
    if not jun_start_dt[chr_name].has_key(leftpos):
        jun_start_dt[chr_name][leftpos] = set()
    jun_start_dt[chr_name][leftpos].add(rightpos)

    if not jun_end_dt.has_key(chr_name):
        jun_end_dt[chr_name] = {}
    if not jun_end_dt[chr_name].has_key(rightpos):
        jun_end_dt[chr_name][rightpos] = set()
    jun_end_dt[chr_name][rightpos].add(leftpos)

jun_file.close()

####Add junctions from RefSeq###################################################

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

    if not jun_start_dt.has_key(chr_name):
        jun_start_dt[chr_name] = {}
    if not jun_end_dt.has_key(chr_name):
        jun_end_dt[chr_name] = {}

    for jun_start in Exon_end_list[:-1]:
        jun_start = int(jun_start)
        jun_end = int(Exon_start_list[i])

        if not jun_start_dt[chr_name].has_key(jun_start):
             jun_start_dt[chr_name][jun_start] = set()
        jun_start_dt[chr_name][jun_start].add(jun_end)

        if not jun_end_dt[chr_name].has_key(jun_end):
            jun_end_dt[chr_name][jun_end] = set()
        jun_end_dt[chr_name][jun_end].add(jun_start)

        locus = chr_name + ':' + str(jun_start) + '-' +str(jun_end)
        jun_strand_dt[locus] = strand

        i += 1

for chr_name in jun_start_dt:
    for jun_start in jun_start_dt[chr_name]:
        jun_start_dt[chr_name][jun_start] = list(jun_start_dt[chr_name][jun_start])      
for chr_name in jun_end_dt:
    for jun_end in jun_end_dt[chr_name]:
        jun_end_dt[chr_name][jun_end] = list(jun_end_dt[chr_name][jun_end])

#get jun for gene###############################################################

exon_boundary_dt = {}
gene_region_dt = {}
exon_boundary_file = open(exon_boundary_filename,'r')
for line in exon_boundary_file:
    ls = line.strip().split('\t')
    locus = ls[0]
    locus_ls = locus.split(':')
    locus_ls_ls = locus_ls[1].split('-')
    start = int(locus_ls_ls[0])
    end = int(locus_ls_ls[1])
    chr_name = locus_ls[0]

    exon_boundary_str_ls = ls[1:]
    exon_boundary_ls = []
    for exon_boundary in exon_boundary_str_ls:
        exon_boundary_ls.append( int(exon_boundary) )
    exon_boundary_ls.sort()

    exon_boundary_dt[locus]= exon_boundary_ls

    if not gene_region_dt.has_key(chr_name):
        gene_region_dt[chr_name] = {}
    if not gene_region_dt[chr_name].has_key(start):
        gene_region_dt[chr_name][start] = []

    gene_region_dt[chr_name][start].append(end)

exon_boundary_file.close()

#####

def nonredundant_left_threeend(exon_boundary_ls,threeend_ls,gene_start,gene_end):
    result = []
    threeend_ls.sort()
    gene_start_index = searchsorted(threeend_ls,gene_start)
    leftmost_pt_ls = threeend_ls[:gene_start_index]

#    if min(threeend_index_s) == 0: use left most

    if len(leftmost_pt_ls)>0:
        leftmost_pt = min(leftmost_pt_ls)
        if gene_start - leftmost_pt <=10:
            leftmost_pt = gene_start
        result.append( leftmost_pt )

    threeend_index_ls = searchsorted(exon_boundary_ls,threeend_ls,side='right')
    threeend_index_s = set(threeend_index_ls)

    for a in threeend_index_s:
        if a > 0:
            result.append( exon_boundary_ls[a-1] )

    result_dt = {}
    k = 0
    for a in threeend_index_ls:
        pt = threeend_ls[k]
        k +=1
        if a == 0:
            result_dt[pt] = leftmost_pt
        else:
            result_dt[pt] = exon_boundary_ls[a-1] 

    return result_dt, list(set(result))

def nonredundant_right_threeend(exon_boundary_ls,threeend_ls,gene_start,gene_end):
    result = []
    threeend_ls.sort()
    gene_end_index = searchsorted(threeend_ls,gene_end,side='right')
    rightmost_pt_ls = threeend_ls[gene_end_index:]

#    if min(threeend_index_s) == len(exon_boundary_ls): use right most

    if len(rightmost_pt_ls) >0:
        rightmost_pt = max(rightmost_pt_ls)
        if rightmost_pt - gene_end <= 10:
            rightmost_pt = gene_end
        result.append(rightmost_pt)

    threeend_index_ls = searchsorted(exon_boundary_ls,threeend_ls)
    threeend_index_s = set(threeend_index_ls)

    for a in threeend_index_s:
        if a < len(exon_boundary_ls):
            result.append( exon_boundary_ls[a] )

    result_dt = {}
    k = 0
    for a in threeend_index_ls:
        pt = threeend_ls[k]
        k +=1
        if a == len(exon_boundary_ls):
            result_dt[pt] = rightmost_pt
        else:
            result_dt[pt] = exon_boundary_ls[a] 

    return result_dt, list(set(result))

threeend_pt_dt = {}
LR_gpd = open(nonredundant_LR_gpd_filename ,'r')
for line in LR_gpd:
    ls = line.strip().split('\t')
    chr_name = ls[2]
    gene_name = ls[1]
    strand = ls[3]
    exon_start_list=ls[9].strip(',').split(',')
    exon_end_list=ls[10].strip(',').split(',')
    I_3end, L_3end = check_3end(ls[0])
    if not threeend_pt_dt.has_key(gene_name):
        threeend_pt_dt[gene_name] = [[],[]]
    if I_3end == '-':
        threeend = int(exon_start_list[0])-L_3end
        threeend_pt_dt[gene_name][0].append(threeend)
    elif I_3end == '+':
        threeend = int(exon_end_list[-1])+L_3end
        threeend_pt_dt[gene_name][1].append(threeend)
LR_gpd.close()

ref_threeend_dt = {}
ref_fiveend_dt = {}

if I_ref3end_isoformconstruction == "1":
    ref = open(refFlat_filename,'r')
    for line in ref:
        ls = line.strip().split('\t')
        chr_name = ls[2]
        if not ref_threeend_dt.has_key(chr_name):
            ref_threeend_dt[chr_name] = {}
        if ls[3] == "+":
            ref_threeend_dt[chr_name][ int(ls[5]) ] = ls[3]
        elif ls[3] == "-":
            ref_threeend_dt[chr_name][ int(ls[4]) ] = ls[3]
    ref.close()

if I_ref5end_isoformconstruction == "1":
    ref = open(refFlat_filename,'r')
    for line in ref:
        ls = line.strip().split('\t')
        chr_name = ls[2]
        if not ref_fiveend_dt.has_key(chr_name):
            ref_fiveend_dt[chr_name] = {}        
        if ls[3] == "+":
            ref_fiveend_dt[chr_name][ int(ls[4]) ] = ls[3]
        elif ls[3] == "-":
            ref_fiveend_dt[chr_name][ int(ls[5]) ] = ls[3]
    ref.close()


for chr_name in gene_region_dt:
    if not ref_threeend_dt.has_key(chr_name):
        continue
    gene_region_startls = gene_region_dt[chr_name].keys()
    gene_region_startls.sort()
    ref_threeend_ls = ref_threeend_dt[chr_name].keys()
    ref_threeend_indexls = searchsorted(gene_region_startls,ref_threeend_ls,side='right') 
    i = 0
    for ref_threeend_index in ref_threeend_indexls:
        ref_threeend = ref_threeend_ls[i]
        i += 1
        if ref_threeend_index == 0:
            continue
        gene_region_start = gene_region_startls[ref_threeend_index - 1]
        gene_region_end_ls = gene_region_dt[chr_name][gene_region_start]

        for gene_region_end in gene_region_end_ls:

            if ref_threeend > gene_region_end:
                continue
            locus = chr_name + ":" + str(gene_region_start) + "-" + str(gene_region_end)

            if not threeend_pt_dt.has_key(locus):
                threeend_pt_dt[locus] = [[],[]]
            if ref_threeend_dt[chr_name][ref_threeend] == "-":
                threeend_pt_dt[locus][0].append(ref_threeend)
            elif ref_threeend_dt[chr_name][ref_threeend] == "+":
                threeend_pt_dt[locus][1].append(ref_threeend)

threeend_nonredpt_dt={}
threeend_boundary_dt={}
for gene_name in threeend_pt_dt:
    gene_name_ls=gene_name.split(':')
    gene_name_ls1_ls = gene_name_ls[1].split('-')
    gene_start = int(gene_name_ls1_ls[0])
    gene_end = int(gene_name_ls1_ls[1])

    threeend_nonredpt_dt[gene_name] = [[],[]]
    threeend_boundary_dt[gene_name] = [[],[]]
    exon_boundary_ls = exon_boundary_dt[gene_name]

    threeend_boundary_dt[gene_name][0], threeend_nonredpt_dt[gene_name][0] = nonredundant_left_threeend(exon_boundary_ls,threeend_pt_dt[gene_name][0],gene_start,gene_end)
    threeend_boundary_dt[gene_name][1], threeend_nonredpt_dt[gene_name][1] = nonredundant_right_threeend(exon_boundary_ls,threeend_pt_dt[gene_name][1],gene_start,gene_end)

example_gene = "chr19:21541734-21571384"

#print threeend_pt_dt[example_gene]
#print threeend_nonredpt_dt[example_gene]
#print threeend_boundary_dt[example_gene]

for gene_name in threeend_nonredpt_dt:
    gene_name_ls=gene_name.split(':')
    chr_name = gene_name_ls[0]
    temp_ls = []
    for pt in threeend_nonredpt_dt[gene_name][0]:
        if not jun_start_dt[chr_name].has_key(pt):
            temp_ls.append(pt)
    threeend_nonredpt_dt[gene_name][0] = temp_ls

    temp_ls = []
    for pt in threeend_nonredpt_dt[gene_name][1]:
        if not jun_end_dt[chr_name].has_key(pt):
            temp_ls.append(pt)
    threeend_nonredpt_dt[gene_name][1] = temp_ls    
    

#print threeend_pt_dt[example_gene]
#print threeend_nonredpt_dt[example_gene]
#print threeend_boundary_dt[example_gene]


cap_dt = parse_5cap_file(CAGE_filename)
fwd_5end_dt,rev_5end_dt = define_5end(cap_dt,1)
print "yue", cap_dt
print "kinfai",  fwd_5end_dt,rev_5end_dt
for chr_name in gene_region_dt:
    if not ref_fiveend_dt.has_key(chr_name):
        continue
    gene_region_startls = gene_region_dt[chr_name].keys()
    gene_region_startls.sort()
    ref_fiveend_ls = ref_fiveend_dt[chr_name].keys()
    ref_fiveend_indexls = searchsorted(gene_region_startls,ref_fiveend_ls,side='right')
    i = 0
    for ref_fiveend_index in ref_fiveend_indexls:
        ref_fiveend = ref_fiveend_ls[i]
        i += 1
        if ref_fiveend_index == 0:
            continue
        gene_region_start = gene_region_startls[ref_fiveend_index - 1]
        gene_region_end_ls = gene_region_dt[chr_name][gene_region_start]

        for gene_region_end in gene_region_end_ls:
            if ref_fiveend > gene_region_end:
                continue
            locus = chr_name + ":" + str(gene_region_start) + "-" + str(gene_region_end)
            if not rev_5end_dt.has_key(locus):
                rev_5end_dt[locus]=[]
            if not fwd_5end_dt.has_key(locus):
                fwd_5end_dt[locus]=[]

            if ref_fiveend_dt[chr_name][ref_fiveend] == "-":
                rev_5end_dt[locus].append(ref_fiveend)
            elif ref_fiveend_dt[chr_name][ref_fiveend] == "+":
                fwd_5end_dt[locus].append(ref_fiveend)

for locus in rev_5end_dt:
    rev_5end_dt[locus] = list( set(rev_5end_dt[locus]) )

for locus in fwd_5end_dt:
    fwd_5end_dt[locus] = list( set(fwd_5end_dt[locus]) )    

#print "5end"
#print rev_5end_dt[example_gene]
#print fwd_5end_dt[example_gene]

######

gene_jun_dt = {}
for chr_name in gene_region_dt:
    if not jun_start_dt.has_key(chr_name):
        continue
    gene_start_ls = gene_region_dt[chr_name].keys()
    gene_start_ls.sort()
    jun_start_ls = jun_start_dt[chr_name].keys()
    jun_start_ls.sort() 
    jun_start_index_ls = searchsorted(gene_start_ls,jun_start_ls,side='right')
    i = 0
    for jun_start in jun_start_ls:
        jun_start_index = jun_start_index_ls[i]
        i +=1
        if jun_start_index == 0:
            continue
        gene_start = gene_start_ls[  jun_start_index -1 ]
        if len( gene_region_dt[chr_name][gene_start] ) > 1:
            print "warning", "multi gene end"
            continue
        gene_end = gene_region_dt[chr_name][gene_start][0]
        locus = chr_name + ":" + str(gene_start) + "-" + str(gene_end)

        for jun_end in jun_start_dt[chr_name][jun_start]:
            if gene_end < jun_end:
#kinfai#                print "junction out 3 end"
                continue
            if not gene_jun_dt.has_key(locus):
                gene_jun_dt[locus] = []
            gene_jun_dt[locus].append([jun_start,jun_end])

gene_jun_end_dt = {}
for locus in gene_jun_dt:
    temp_dt = {}
    gene_jun_end_dt[locus] = []
    for jun in gene_jun_dt[locus]:
        if not temp_dt.has_key(jun[1]):
            temp_dt[jun[1]] = []
        temp_dt[jun[1]].append(jun[0])
    jun_end_ls = temp_dt.keys()
    jun_end_ls.sort()
    jun_end_ls.reverse()
    for jun_end in jun_end_ls:
        for jun_start in temp_dt[jun_end]:
            gene_jun_end_dt[locus].append([jun_start,jun_end])
    
################################################################################
def printjunlk(item):
    temp_ls = []
    if len(item) == 1 and item[0] == []:
        return ''
    elif len(item) == 0:
        return ''
    for jun in item:
        temp_ls.append(str(jun[0])+'_'+str(jun[1]))
    return '_'.join(temp_ls)

def filter_strand(gpd_strand, jun_strand_dt, chr_name, jun_ls_ls):
    result = []
    for jun_ls in jun_ls_ls:
        jun_locus = chr_name + ':' + str(jun_ls[0]) + '-' +str(jun_ls[1])
        jun_strand = ''
        if jun_strand_dt.has_key (jun_locus):
            jun_strand = jun_strand_dt[jun_locus]
        if jun_strand == gpd_strand :
            result.append(jun_ls)
    return result
################################################################################

def make_jun_combination(ref_start_lk_ls, target,exon_len_limit):

#    ref_start_lk_ls = [[]]

    for pt in target:
        temp_start_lk_ls = []
        I = 0

        if ref_start_lk_ls == [[]] :
            I = 1
        elif ref_start_lk_ls[-1][-1][1] < pt[0]:
            I = 1

        if I ==0:
           ref_start_lk_ls.append([pt])
        elif I == 1:
            temp_start_lk_ls.extend(ref_start_lk_ls)
            ref_start_lk_ls = deepcopy(ref_start_lk_ls)

            for item in ref_start_lk_ls:

                if item == [] :
                    item.append(pt)
                    temp_start_lk_ls.append(item)
                elif  pt[0] < item[-1][1] + exon_len_limit:
                    item.append(pt)
                    temp_start_lk_ls.append(item)

            ref_start_lk_ls = temp_start_lk_ls

    return ref_start_lk_ls

def make_endjun_combination(ref_start_lk_ls, target,exon_len_limit):
#    ref_start_lk_ls = [[[1,6]]]
    for pt in target:
        temp_start_lk_ls = []
        temp_start_lk_ls.extend(ref_start_lk_ls)
        ref_start_lk_ls = deepcopy(ref_start_lk_ls)
        for item in ref_start_lk_ls:
            if item[-1][1] < pt[0] and pt[0] < item[-1][1] + exon_len_limit:
                item.append(pt)
                temp_start_lk_ls.append(item)
        ref_start_lk_ls = temp_start_lk_ls
    for item in ref_start_lk_ls:
        del item[0]
    return ref_start_lk_ls

def make_startjun_combination(ref_start_lk_ls, target,exon_len_limit):
#    ref_start_lk_ls = [[[1,6]]]
#    target.reverse()
    for pt in target:
        temp_start_lk_ls = []
        temp_start_lk_ls.extend(ref_start_lk_ls)
        ref_start_lk_ls = deepcopy(ref_start_lk_ls)
        for item in ref_start_lk_ls:
            if item[-1][0] > pt[1] and pt[1] > item[-1][0] - exon_len_limit:
                item.append(pt)
                temp_start_lk_ls.append(item)
        ref_start_lk_ls = temp_start_lk_ls
    for item in ref_start_lk_ls:
        del item[0]
        item.reverse()
    return ref_start_lk_ls

####

def make_endjun_threeend_combination(old_ref_start_lk_ls, target,threeend_nonredpt_ls,exon_len_limit,exon_boundary_ls, raw_ending_pt):
#    ref_start_lk_ls = [[[1,6]]]
    result = {}
    for threeend_pt in threeend_nonredpt_ls:
        ref_start_lk_ls = deepcopy(old_ref_start_lk_ls)
        result[threeend_pt] = []
        temp_target_ls = []
        temp_result = []

        if ref_start_lk_ls[-1][-1][1] > threeend_pt:
            continue

        for pt in target:
            if ref_start_lk_ls[-1][-1][1]<pt[0] and pt[1] < threeend_pt:  ###
                temp_target_ls.append(pt)
        for pt in temp_target_ls:
            temp_start_lk_ls = []
            temp_start_lk_ls.extend(ref_start_lk_ls)
            ref_start_lk_ls = deepcopy(ref_start_lk_ls)
            for item in ref_start_lk_ls:
                if item[-1][1] < pt[0] and pt[0] < item[-1][1] + exon_len_limit:  ####
                    item.append(pt)
                    temp_start_lk_ls.append(item)
                    if threeend_pt - pt[1] < exon_len_limit: ###
                        temp_result.append(item)
            ref_start_lk_ls = temp_start_lk_ls

#        print temp_target_ls
        if temp_target_ls == []:
            temp_target_ls = old_ref_start_lk_ls
#            print temp_target_ls, "yyy"

        for item in ref_start_lk_ls:
            if item[-1][1] > threeend_pt - exon_len_limit: ###
                del item[0]
                result[threeend_pt].append(item)

        temp_ending_index_ls = searchsorted(exon_boundary_ls, [raw_ending_pt,threeend_pt])
        if temp_ending_index_ls[0] == temp_ending_index_ls[1]:
            result[threeend_pt].append([])

    result_ls = []
    ending_ls = []
    for threeend_pt in result:
        for item in result[threeend_pt]:
            ending_ls.append(threeend_pt)
            result_ls.append(item)
#    print ending_ls
#    print result_ls

    return     ending_ls, result_ls

def make_startjun_threeend_combination(old_ref_start_lk_ls, target,threeend_nonredpt_ls,exon_len_limit,exon_boundary_ls, raw_ending_pt):
#    target.reverse()
    
    result = {}
    for threeend_pt in threeend_nonredpt_ls:
        ref_start_lk_ls = deepcopy(old_ref_start_lk_ls)
        result[threeend_pt] = []
        temp_target_ls = []
        temp_result = []

        if ref_start_lk_ls[0][0][0] < threeend_pt:
            continue

        for pt in target:
            if ref_start_lk_ls[0][0][0] > pt[1] and pt[0] > threeend_pt:
                temp_target_ls.append(pt)
        for pt in temp_target_ls:
            temp_start_lk_ls = []
            temp_start_lk_ls.extend(ref_start_lk_ls)
            ref_start_lk_ls = deepcopy(ref_start_lk_ls)
            for item in ref_start_lk_ls:
                if item[-1][0] > pt[1] and pt[1] > item[-1][0] - exon_len_limit:
                    item.append(pt)
                    temp_start_lk_ls.append(item)
                    if pt[0] - threeend_pt < exon_len_limit: ###
                        temp_result.append(item)
            ref_start_lk_ls = temp_start_lk_ls

        if temp_target_ls == []:
            temp_target_ls = old_ref_start_lk_ls

        for item in ref_start_lk_ls:
            if item[-1][0] - threeend_pt < exon_len_limit: ###
                del item[0]
                item.reverse()
                result[threeend_pt].append(item)

        temp_ending_index_ls = searchsorted(exon_boundary_ls, [raw_ending_pt,threeend_pt])
        if temp_ending_index_ls[0] == temp_ending_index_ls[1]:
            result[threeend_pt].append([])

    result_ls = []
    ending_ls = []

    for threeend_pt in result:
        for item in result[threeend_pt]:
            ending_ls.append(threeend_pt)
            result_ls.append(item)
    return     ending_ls, result_ls

####

LR_gpd_dt = {}
candidate_dt = {}
LR_gpd = open(nonredundant_LR_gpd_filename ,'r')

i_LR_gpd = 1
for line in LR_gpd:
    ls = line.strip().split('\t')
    chr_name = ls[2]
    locus = ls[1]
    exon_boundary_ls = exon_boundary_dt[locus]
    exon_start_list=ls[9].strip(',').split(',')
    exon_end_list=ls[10].strip(',').split(',')

    LR_start = int(exon_start_list[0])
    LR_end = int(exon_end_list[-1])

    i = 0
    jun_ls = []
    gpd_strand_ls =['']
    for jun_end in exon_start_list[1:]:
        jun_start = exon_end_list[i]

        jun_locus = chr_name + ':' + jun_start + '-' + jun_end
        if jun_strand_dt.has_key(jun_locus):
            gpd_strand_ls.append( jun_strand_dt[jun_locus] )

        jun_ls.append(jun_start)
        jun_ls.append(jun_end)
        i += 1
    jun_lk = '_'.join(jun_ls)

    gpd_strand_s = set(gpd_strand_ls)
    L_gpd_strand_s = len(gpd_strand_s)
    if L_gpd_strand_s == 1:
        gpd_strand = ''
    elif L_gpd_strand_s == 2:
        gpd_strand = list( gpd_strand_s - set(['']) )[0] 
    else:
        gpd_strand = ''
        print "gene fusion:", line.strip()

    first_jun = [int(exon_end_list[0]), int(exon_start_list[1]) ]
    last_jun = [int(exon_end_list[-2]), int(exon_start_list[-1]) ]

    chr_name = ls[2]
    locus = ls[1]

    locus_ls = locus.split(":")
    chr_name = locus_ls[0]

    locus_ls1 = locus_ls[1].split("-")
    gene_start = locus_ls1[0] 
    gene_end = locus_ls1[1]

#    if not LR_gpd_dt.has_key(chr_name):
#        LR_gpd_dt[chr_name]={}
#    if not LR_gpd_dt[chr_name].has_key(exon_start_list[0]):
#        LR_gpd_dt[chr_name][exon_start_list[0]]= {}
#    LR_gpd_dt[chr_name][exon_start_list[0]][exon_end_list[-1]] = jun_lk

    start_jun_ls = []
    end_jun_ls = []

    if gene_jun_dt.has_key(locus):
        for item in gene_jun_dt[locus]:
            if item[0] >= LR_end:               #            if item[1]<= LR_start:
                end_jun_ls.append(item)         #                start_jun_ls.append(item)

    if gene_jun_end_dt.has_key(locus):
        for item in gene_jun_end_dt[locus]:
            if item[1]<= LR_start:
                start_jun_ls.append(item)

    start_jun_ls = filter_strand(gpd_strand, jun_strand_dt, chr_name, start_jun_ls)
    end_jun_ls = filter_strand(gpd_strand, jun_strand_dt, chr_name, end_jun_ls)

    I_stop = 0    
    error_str = ""

    if len(start_jun_ls) > limit:
        error_str =  error_str + "len(start_jun_ls) " + str(len(start_jun_ls)) + " > " + str(limit) 
        I_stop = 1
    if len(end_jun_ls) > limit :
        if I_stop == 1:
            error_str =  error_str + " "
        error_str =  error_str + "len(end_jun_ls)" + str(len(end_jun_ls)) + " > " + str(limit)

        I_stop = 1

    if  I_stop == 1 :
        print "warning", error_str + "\t" + locus
        continue

    I_3end, L_3end = check_3end(ls[0])
    if I_3end != gpd_strand:
        I_3end = ''

    ref_start_lk_ls = [[first_jun]]

#    junlk_set = set()

    if not candidate_dt.has_key(locus):
        candidate_dt[locus] = {}
    
    start_junlk_str_ls = []
    end_junlk_str_ls = []

    if I_3end == '+':
        ending_boundary = threeend_boundary_dt[locus][1][LR_end + L_3end]
        
        if not fwd_5end_dt.has_key(locus):
            start_jun_lk_ls = make_startjun_combination( ref_start_lk_ls, start_jun_ls,exon_len_limit)
            left_ending_ls = [-1]*len(start_jun_lk_ls)
        else:
            left_ending_ls, start_jun_lk_ls = make_startjun_threeend_combination(ref_start_lk_ls, start_jun_ls,fwd_5end_dt[locus],exon_len_limit,exon_boundary_ls, exon_start_list[0])

        for item in start_jun_lk_ls:
            start_junlk_str_ls.append( printjunlk(item) )

        p = 0
        for start_junlk_str in start_junlk_str_ls:
            fiveend_ending_str = str(left_ending_ls[p])
            p+=1
            item = start_junlk_str + '_' + jun_lk
            item = item.strip('_')
            if not candidate_dt[locus].has_key(item):
                candidate_dt[locus][item] = {}
            if not candidate_dt[locus][item].has_key(gpd_strand):
                candidate_dt[locus][item][gpd_strand] = set()

            temp_end = int( item.split('_')[-1] ) 
            if ending_boundary - temp_end > exon_len_limit:
                ending_boundary = -1 ####LR_end + L_3end

            candidate_dt[locus][item][gpd_strand].add( fiveend_ending_str + '_' + str(ending_boundary) )
        end_junlk_str_ls = ['']
#            junlk_set.add( start_junlk_str + '-' + jun_lk )
    
    if I_3end == '-':

        ending_boundary = threeend_boundary_dt[locus][0][LR_start - L_3end]

        if not rev_5end_dt.has_key(locus):
            end_jun_lk_ls = make_endjun_combination([[last_jun]],end_jun_ls,exon_len_limit)
            right_ending_ls = [-1]*len(end_jun_lk_ls)
        else:
            right_ending_ls, end_jun_lk_ls = make_endjun_threeend_combination([[last_jun]], end_jun_ls,rev_5end_dt[locus],exon_len_limit,exon_boundary_ls, exon_end_list[-1])

        for item in end_jun_lk_ls:
            end_junlk_str_ls.append( printjunlk(item) )

        p = 0
        for end_junlk_str in end_junlk_str_ls:
            fiveend_ending_str = str(right_ending_ls[p])
            p+=1

            item = jun_lk + '_' + end_junlk_str
            item = item.strip('_')
            if not candidate_dt[locus].has_key(item):
                candidate_dt[locus][item] = {}
            if not candidate_dt[locus][item].has_key(gpd_strand):
                candidate_dt[locus][item][gpd_strand] = set()

            temp_start = int( item.split('_')[0] ) 
            if temp_start - ending_boundary > exon_len_limit:
                ending_boundary = -1 ###LR_start - L_3end

            candidate_dt[locus][item][gpd_strand].add( str(ending_boundary) + '_' + fiveend_ending_str )
        start_junlk_str_ls = ['']
#            junlk_set.add( jun_lk + '-' + end_junlk_str )


    if I_3end == '':

        if threeend_nonredpt_dt.has_key(locus):

            if gpd_strand == '+':
                right_ending_ls, end_jun_lk_ls = make_endjun_threeend_combination([[last_jun]], end_jun_ls,threeend_nonredpt_dt[locus][1],exon_len_limit,exon_boundary_ls, exon_end_list[-1])
                if not fwd_5end_dt.has_key(locus):
                    start_jun_lk_ls = make_startjun_combination( ref_start_lk_ls, start_jun_ls,exon_len_limit)
                    if start_jun_lk_ls ==[]:
                        start_jun_lk_ls = [[]]
                    left_ending_ls = [-1]*len(start_jun_lk_ls)
                else:
                    left_ending_ls, start_jun_lk_ls = make_startjun_threeend_combination(ref_start_lk_ls, start_jun_ls,fwd_5end_dt[locus],exon_len_limit,exon_boundary_ls, exon_start_list[0])
#                    print "KKKK",left_ending_ls, start_jun_lk_ls , right_ending_ls, end_jun_lk_ls 

            elif gpd_strand == '-':
                left_ending_ls, start_jun_lk_ls = make_startjun_threeend_combination(ref_start_lk_ls, start_jun_ls,threeend_nonredpt_dt[locus][0],exon_len_limit,exon_boundary_ls, exon_start_list[0])
                if not rev_5end_dt.has_key(locus):
                    end_jun_lk_ls = make_endjun_combination([[last_jun]],end_jun_ls,exon_len_limit)
                    if end_jun_lk_ls ==[]:
                        end_jun_lk_ls = [[]]
                    right_ending_ls = [-1]*len(end_jun_lk_ls)
                else:
                    right_ending_ls, end_jun_lk_ls = make_endjun_threeend_combination([[last_jun]], end_jun_ls,rev_5end_dt[locus],exon_len_limit,exon_boundary_ls, exon_end_list[-1])

            else:
                start_jun_lk_ls = make_startjun_combination( ref_start_lk_ls, start_jun_ls,exon_len_limit)
                end_jun_lk_ls = make_endjun_combination([[last_jun]],end_jun_ls,exon_len_limit)
                left_ending_ls = [-1]*len(start_jun_lk_ls)
                right_ending_ls = [-1]*len(end_jun_lk_ls)

        else:
            start_jun_lk_ls = make_startjun_combination( ref_start_lk_ls, start_jun_ls,exon_len_limit)
            end_jun_lk_ls = make_endjun_combination([[last_jun]],end_jun_ls,exon_len_limit)
            left_ending_ls = [-1]*len(start_jun_lk_ls)
            right_ending_ls = [-1]*len(end_jun_lk_ls)

        for item in start_jun_lk_ls:
            start_junlk_str_ls.append( printjunlk(item) )
        for item in end_jun_lk_ls:
            end_junlk_str_ls.append( printjunlk(item) )

        m = 0
        for start_junlk_str in start_junlk_str_ls:
            left_ending_boundary = left_ending_ls[m]
            n = 0
            for end_junlk_str in end_junlk_str_ls:
                right_ending_boundary = right_ending_ls[n]

                if gpd_strand != '':
                   item = start_junlk_str + '_' + jun_lk + '_' + end_junlk_str
                else:
                   item = start_junlk_str + '_' + jun_lk + '_' + end_junlk_str
                item = item.strip('_')

                if not candidate_dt[locus].has_key(item):
                    candidate_dt[locus][item] = {}
                if not candidate_dt[locus][item].has_key(gpd_strand):
                    candidate_dt[locus][item][gpd_strand] = set()


                temp_start = int( item.split('_')[0] ) 
                if left_ending_boundary != -1:
                    if temp_start - left_ending_boundary > exon_len_limit:
                        left_ending_boundary = -1

                temp_end = int( item.split('_')[-1] ) 
                if right_ending_boundary != -1:
                    if right_ending_boundary - temp_end > exon_len_limit:
                        right_ending_boundary = -1

                twoending_str = str(left_ending_boundary) + '_' + str(right_ending_boundary)
                candidate_dt[locus][item][gpd_strand].add( twoending_str )
                n+=1
            m+=1
#                    junlk_set.add( start_junlk_str + '-' + jun_lk + '-' + end_junlk_str )

#    if not candidate_dt.has_key(locus):
#        candidate_dt[locus] = set()
#    candidate_dt[locus] = candidate_dt[locus] | junlk_set 

#    if i_LR_gpd%100 == 0:
#        print "jun", i_LR_gpd
#    i_LR_gpd+=1

LR_gpd.close()

###########################################################################################################
output = open(output_filename,'w')
for locus in candidate_dt:
    for item in candidate_dt[locus]:
        for gpd_strand in candidate_dt[locus][item]:
            for twoending_str in candidate_dt[locus][item][gpd_strand]:
                output.write( 'kinfai'.join([locus, item, gpd_strand, twoending_str]) + '\n' )
output.close()
###########################################################################################################



################################################################################

