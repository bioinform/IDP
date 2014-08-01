#!/usr/bin/python


import sys
import os
from numpy import *

if len(sys.argv) >= 3:
    LR_gpd_filename = sys.argv[1]
    exon_boundary_filename = sys.argv[2]
else:
    print("usage: ./compatible_LR.py LR_gpd_filename exon_boundary_filename")
    print("or ")
    sys.exit(1)
################################################################################

def check_jun(jun_list1,jun_list2):
    if (not jun_list1[1]<jun_list2[0]) and (not jun_list2[1]<jun_list1[0]):
        return 1
    return 0
################################################################################

exon_boundary_dt = {}
gene_region_dt = {}
exon_boundary = open(exon_boundary_filename,'r')
for line in exon_boundary:
    ls = line.strip().split('\t')
    locus = ls[0]
    locus_ls = locus.split(':')
    locus_ls_ls = locus_ls[1].split('-')
    start = int(locus_ls_ls[0])
    end = int(locus_ls_ls[1])
    chr_name = locus_ls[0]
    
    if not exon_boundary_dt.has_key(chr_name):
        exon_boundary_dt[chr_name] = {}
    exon_boundary_dt[chr_name][start]=[end,ls[1:]]
    if not gene_region_dt.has_key(chr_name):
        gene_region_dt[chr_name] = []
    gene_region_dt[chr_name].append(start)

exon_boundary.close()

################################################################################

LR_gpd_dict = {}
LR_gpd_start_dict = {}

LR_gpd = open(LR_gpd_filename,'r')
for line in LR_gpd:
    line_list=line.split()
    Exon_start=int(line_list[4])
    Exon_end=int(line_list[5])
    Exon_start_list=line_list[9].strip(",").split(',')
    Exon_end_list=line_list[10].strip(",").split(',')
    strand=line_list[3]
    chr_name = line_list[2]

    start = int(Exon_start_list[0])
    end = int(Exon_end_list[-1])

    if not LR_gpd_start_dict.has_key(chr_name):
        LR_gpd_start_dict[chr_name]=[]
    LR_gpd_start_dict[chr_name].append(start)

    if not LR_gpd_dict.has_key(chr_name):
        LR_gpd_dict[chr_name]=[]
    LR_gpd_dict[chr_name].append([end,line_list,int(line_list[8])])

LR_gpd.close()

gene_LR_gpd_dt={}
for chr_name in LR_gpd_start_dict:
    if not gene_region_dt.has_key(chr_name):
        continue
    gene_LR_gpd_dt[chr_name] = {}
    gene_region_dt[chr_name].sort()
    LR_gpd_index_ls = searchsorted(gene_region_dt[chr_name],LR_gpd_start_dict[chr_name],side='right')
    N_gene = len(gene_region_dt[chr_name])

    i  = 0

    for LR_gpd_index in LR_gpd_index_ls:

        temp_LR_gpd = LR_gpd_dict[chr_name][i]
        temp_start = LR_gpd_start_dict[chr_name][i]

        i += 1

        gene_start = gene_region_dt[chr_name][LR_gpd_index-1]
        gene_end = exon_boundary_dt[chr_name][gene_start][0]

        if LR_gpd_index == 0 :
            if gene_start > temp_LR_gpd[0]:
                print "5' terminus chrom", temp_LR_gpd
                continue      ########################
            elif gene_start <= temp_LR_gpd[0]:
                print "5' terminus chrom intersected with gene", temp_LR_gpd

        if gene_end < temp_LR_gpd[0]:

            if N_gene == LR_gpd_index :
                print "terminus chrom:\t", temp_LR_gpd
                continue      ########################
        
            else:
                next_gene_start = gene_region_dt[chr_name][LR_gpd_index]
                next_gene_end = exon_boundary_dt[chr_name][next_gene_start][0]

                if check_jun(  [next_gene_start,next_gene_end], [temp_start,temp_LR_gpd[0]]  ) == 1:
                    if not gene_LR_gpd_dt[chr_name].has_key(next_gene_start):
                        gene_LR_gpd_dt[chr_name][next_gene_start]=[]
                    gene_LR_gpd_dt[chr_name][next_gene_start].append(temp_LR_gpd)
                    continue         ########################
                print "3'end gene:\t", temp_LR_gpd
                continue

        if not gene_LR_gpd_dt[chr_name].has_key(gene_start):
            gene_LR_gpd_dt[chr_name][gene_start]=[]
        gene_LR_gpd_dt[chr_name][gene_start].append(temp_LR_gpd)


#print "kinfai", gene_LR_gpd_dt["chr8"][74702841]

################################################################################
def make_jun_str(Exon_start_list_str,Exon_end_list_str):
    Exon_start_list = Exon_start_list_str.strip(',').split(',')
    Exon_end_list = Exon_end_list_str.strip(',').split(',')

    jun_ls = []
    i = 0
    for jun_end in Exon_start_list[1:]: 
        jun_start = Exon_end_list[i]
        jun_ls.append(jun_start)
        jun_ls.append(jun_end)
        i += 1
    return '_'.join(jun_ls)

def comare_ls(ls1,ls2):
#    return 0 if not identical
#    return 1 if identical

    i = 0
    for element in ls1:
        if element != ls2[i]:
            return 0
        i += 1
    return 1

def calculatestat(item):
    I = item.find('-')
    if I == -1:
        I = item.find('+')
        if I == -1:
            stat_str = item
        else:
            stat_str = item[:I]
    else:
        stat_str = item[:I]

    item_ls = stat_str.split("|")[0].split("_")
    return float(item_ls[0])*float(item_ls[1])

def check_3end(readname):
#	- is left
#	+ is right
#	"" is no polyA3end
    I = readname.find('-')
    if I == -1:
        I = readname.find('+')
        if I == -1:
            return "",0
        else:
            return '+',int(readname[I+1:])
    else:
        return '-',int(readname[I+1:])


def compatiable(gpd_ls, ref_LR_gpd,exon_boundary_ls):
#    return 0 if compatiable
#    return 1 if incompatiable and ref_LR_gpd is larger
#    return -1 if incompatiable and gpd_ls is larger

    jun_str = gpd_ls[11]
    ref_jun_str = ref_LR_gpd[11]

    if ref_jun_str.find(jun_str) == -1:
        return 0
    else:

        
        I_ref_3end,L_ref_3end = check_3end(ref_LR_gpd[0])
        I_tag_3end,L_tag_3end = check_3end(gpd_ls[0])

        first_exon_start = int(gpd_ls[9].strip(',').split(',')[0])
        last_exon_end = int(gpd_ls[10].strip(',').split(',')[-1])

        ref_first_exon_start = int(ref_LR_gpd[9].strip(',').split(',')[0])
        ref_last_exon_end = int(ref_LR_gpd[10].strip(',').split(',')[-1])


        if I_ref_3end == '-':
            ref_first_exon_start -= L_ref_3end
        elif I_ref_3end == '+':
            ref_last_exon_end += L_ref_3end

        if I_tag_3end == '-':
            first_exon_start -= L_tag_3end
        elif I_tag_3end == '+':
            last_exon_end += L_tag_3end

        terminal_splice_index = searchsorted(exon_boundary_ls,[first_exon_start, last_exon_end, ref_first_exon_start, ref_last_exon_end])

        if terminal_splice_index[0] < terminal_splice_index[2] and terminal_splice_index[1] > terminal_splice_index[3]:
            if I_ref_3end != '' and I_tag_3end == '':
                return 0
            elif I_ref_3end == '' and I_tag_3end != '':
                return -1
            elif I_ref_3end != '' and I_tag_3end != '':
                return 0
            elif I_ref_3end == '' and I_tag_3end == '':
                return -1

        elif terminal_splice_index[0] < terminal_splice_index[2] and terminal_splice_index[1] == terminal_splice_index[3]:
            if I_ref_3end != '' and I_tag_3end == '':
                return 0
            elif I_ref_3end == '' and I_tag_3end != '':
                return -1
            elif I_ref_3end != '' and I_tag_3end != '':
                if I_ref_3end == '+' and I_tag_3end == '+':
                    return -1
                else:
                    return 0
            elif I_ref_3end == '' and I_tag_3end == '':
                return -1

        elif terminal_splice_index[0] < terminal_splice_index[2] and terminal_splice_index[1] < terminal_splice_index[3]:
            return 0

        elif terminal_splice_index[0] == terminal_splice_index[2] and terminal_splice_index[1] > terminal_splice_index[3]:
            if I_ref_3end != '' and I_tag_3end == '':
                return 0
            elif I_ref_3end == '' and I_tag_3end != '':
                return -1
            elif I_ref_3end != '' and I_tag_3end != '':
                if I_ref_3end == '-' and I_tag_3end == '-':
                    return -1
                else:
                    return 0
            elif I_ref_3end == '' and I_tag_3end == '':
                return -1

        elif terminal_splice_index[0] == terminal_splice_index[2] and terminal_splice_index[1] == terminal_splice_index[3]:
            if I_ref_3end != '' and I_tag_3end == '':
                return 1
            elif I_ref_3end == '' and I_tag_3end != '':
                return -1
            elif I_ref_3end != '' and I_tag_3end != '':
                if I_ref_3end == I_tag_3end:
                    if calculatestat(gpd_ls[0]) > calculatestat(ref_LR_gpd[0]):
                        return -1
                    else:
                        return 1
                else:
                    return 0
            elif I_ref_3end == '' and I_tag_3end == '':
                if calculatestat(gpd_ls[0]) > calculatestat(ref_LR_gpd[0]):
                    return -1
                else:
                    return 1

        elif terminal_splice_index[0] == terminal_splice_index[2] and terminal_splice_index[1] < terminal_splice_index[3]:
            if I_ref_3end != '' and I_tag_3end == '':
                return 1
            elif I_ref_3end == '' and I_tag_3end != '':
                return 0
            elif I_ref_3end != '' and I_tag_3end != '':
                if I_ref_3end == '-' and I_tag_3end == '-':
                    return 1
                else:
                    return 0
            elif I_ref_3end == '' and I_tag_3end == '':
                return 1

        elif terminal_splice_index[0] > terminal_splice_index[2] and terminal_splice_index[1] > terminal_splice_index[3]:
            return 0

        elif terminal_splice_index[0] > terminal_splice_index[2] and terminal_splice_index[1] == terminal_splice_index[3]:
            if I_ref_3end != '' and I_tag_3end == '':
                return 1
            elif I_ref_3end == '' and I_tag_3end != '':
                return 0
            elif I_ref_3end != '' and I_tag_3end != '':
                if I_ref_3end == '+' and I_tag_3end == '+':
                    return 1
                else:
                    return 0
            elif I_ref_3end == '' and I_tag_3end == '':
                return 1

        elif terminal_splice_index[0] > terminal_splice_index[2] and terminal_splice_index[1] < terminal_splice_index[3]:
            if I_ref_3end != '' and I_tag_3end == '':
                return 1
            elif I_ref_3end == '' and I_tag_3end != '':
                return 0
            elif I_ref_3end != '' and I_tag_3end != '':
                return 0
            elif I_ref_3end == '' and I_tag_3end == '':
                return 1

def nonredudant(ref_LR_gpd_ls,gpd_ls_ls,exon_boundary_ls):
    for gpd_ls in gpd_ls_ls:

        if ref_LR_gpd_ls == []:
            ref_LR_gpd_ls = [gpd_ls]
#            print ref_LR_gpd_ls

        j = 0
        del_index_ls = []

        for ref_LR_gpd in ref_LR_gpd_ls:
            Icompatiable = compatiable(gpd_ls, ref_LR_gpd,exon_boundary_ls)
##            print gpd_ls
##            print ref_LR_gpd
##            print Icompatiable

            if Icompatiable == 1:
                break
            elif Icompatiable == -1:
                del_index_ls.append(j)
            j += 1    
#        print         del_index_ls
        del_index_ls.reverse()
        for j in del_index_ls:
            del ref_LR_gpd_ls[j]
        if Icompatiable != 1:
            ref_LR_gpd_ls.append(gpd_ls)

    return ref_LR_gpd_ls

def printout(ref_LR_gpd_ls,gene_name):
    if 1:
        for ref_LR_gpd in ref_LR_gpd_ls:
            ref_LR_gpd[1] = gene_name
            print "\t".join(ref_LR_gpd)

################################################################################

def GetFullJun(exon_boundary_ls,item1):
    jun_ls = []
    exon_start_ls = item1[9].strip(',').split(',')
    exon_end_ls = item1[10].strip(',').split(',')
    m = 0
    for exon_start in exon_start_ls:
        exon_end = exon_end_ls[m]
        start_i = searchsorted(exon_boundary_ls,int(exon_start),side='right')
        end_i = searchsorted(exon_boundary_ls,int(exon_end) )

        if len(exon_boundary_ls) > 0 and start_i < len(exon_boundary_ls):
            if m == 0:
                if exon_boundary_ls[start_i] - int(exon_start) <= 10:
                    exon_start_ls[0] = str( exon_boundary_ls[start_i] )
                    start_i += 1

            if m == len(exon_end_ls)-1:
                if int(exon_end) - exon_boundary_ls[end_i-1] <= 10:
                    exon_end_ls[-1] = str( exon_boundary_ls[end_i-1] ) 
                    end_i -= 1

        for pt in exon_boundary_ls[start_i:end_i]:
            jun_ls.append( str(pt)+"_"+str(pt+1) )

        if m < len(exon_end_ls)-1:
            jun_ls.append(  exon_end+"_"+exon_start_ls[m+1]  )
        m += 1

    item1[9] = ','.join(exon_start_ls) + ','
    item1[10] = ','.join(exon_end_ls) + ','
    item1.insert(11,'_'.join(jun_ls) )
    return item1

for chr_name in gene_LR_gpd_dt:
    for gene_start in gene_LR_gpd_dt[chr_name]:
        gene_end = exon_boundary_dt[chr_name][gene_start][0]
        exon_boundary_str_ls = exon_boundary_dt[chr_name][gene_start][1]

        exon_boundary_ls = []
        for exon_boundary in exon_boundary_str_ls:
            exon_boundary_ls.append(int(exon_boundary))

        gene_name = chr_name + ":" + str(gene_start) + "-" +str(gene_end)

        temp_LR_gpd_dt = {}
        for item in gene_LR_gpd_dt[chr_name][gene_start]:

            if not temp_LR_gpd_dt.has_key(item[2]):
                temp_LR_gpd_dt[item[2]]=[]

            temp_item1 = GetFullJun(exon_boundary_ls,item[1])


            temp_LR_gpd_dt[item[2]].append(temp_item1)

        Njun_list = temp_LR_gpd_dt.keys()
        Njun_list.sort()
        Njun_list.reverse()

##        print "################################################################################"
        ref_LR_gpd_ls = []
        for Njun in Njun_list:
##            printout(temp_LR_gpd_dt[Njun])
##            print str(Njun)*20
            ref_LR_gpd_ls = nonredudant(ref_LR_gpd_ls,temp_LR_gpd_dt[Njun],exon_boundary_ls)
##            printout(ref_LR_gpd_ls)


        printout(ref_LR_gpd_ls,gene_name)

