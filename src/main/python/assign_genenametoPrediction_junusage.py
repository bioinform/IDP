#!/usr/bin/python

import sys
import os
from numpy import *
################################################################################

if len(sys.argv) >= 4:
    refFlat_filename = sys.argv[1]
    gpd_filename =  sys.argv[2]
    output_filename =  sys.argv[3]
else:
    print("usage: python refFlat.txt gpd_filename output_filename")
    print("or ./g refFlat.txt gpd_filename output_filename")
    sys.exit(1)

################################################################################

def addrefFlat(ref_iso_dt,ref_dt, ref_refFlat_filename):
    ref=open(ref_refFlat_filename,'r')
    for refline in ref:
        ls = refline.strip().split('\t')
        gene = ls[0]
        ID = ls[1]
        chr_name = ls[2]
        if ls[8]=="1":
            continue
        jun_end_ls = ls[9].strip(',').split(',')[1:]
        jun_start_ls = ls[10].strip(',').split(',')[:-1]
        locus = chr_name+':'+'_'.join(jun_start_ls)+'-'+('_').join(jun_end_ls)
        if not ref_iso_dt.has_key(locus):
            ref_iso_dt[locus]=[[],[]]
        ref_iso_dt[locus][0].append(gene)
        ref_iso_dt[locus][1].append(ID)
        if not ref_dt.has_key(locus):
            ref_dt[locus]=[]
        ref_dt[locus].append(refline)
    ref.close()

################################################################################

ref_iso_dt = {}
ref_dt = {}
addrefFlat(ref_iso_dt,ref_dt, refFlat_filename)
################################################################################

name_dt={}

tag = open(gpd_filename,'r')
for line in tag:
    ls = line.strip().split('\t')
    gene = ls[0]
    ID = ls[1]
    chr_name = ls[2]
    if ls[8]=="1":
        continue
    jun_end_ls = ls[9].strip(',').split(',')[1:]
    jun_start_ls = ls[10].strip(',').split(',')[:-1]
    locus = chr_name+':'+'_'.join(jun_start_ls)+'-'+('_').join(jun_end_ls)

    if ref_iso_dt.has_key(locus):
        if len(ref_iso_dt[locus][0])>1:
            print "mutligene for:", ls[0:2],ref_iso_dt[locus][0]
        name_dt[ID]=[gene,ref_iso_dt[locus][0][0],ref_iso_dt[locus][1]]
        
tag.close()

################################################################################

ref=open(refFlat_filename,'r')

jun_transcript_ID_dict = {}

for refline in ref:

    refline_list=refline.strip().split('\t')

    exon_start_list=refline_list[9].strip(',').split(',')
    exon_end_list=refline_list[10].strip(',').split(',')

    jun_start = array(exon_end_list[:-1],dtype=int)
    jun_end = array(exon_start_list[1:],dtype=int)

    gene_id=refline_list[0]
    transcript_id = refline_list[1]

    chr_name = refline_list[2]

    i = 0
    for jun_end in exon_start_list[1:]:
        jun_start = exon_end_list[i]
        jun = chr_name + ":" + jun_start +  "_" + jun_end

        if not jun_transcript_ID_dict.has_key(jun):
             jun_transcript_ID_dict[jun] = []
        jun_transcript_ID_dict[jun].append(gene_id)
        i+=1

ref.close()

################################################################################

target=open(gpd_filename,'r')

for line in target:
    ls = line.strip().split('\t')

    exon_start_list=ls[9].strip(',').split(',')
    exon_end_list=ls[10].strip(',').split(',')

    jun_start = array(exon_end_list[:-1],dtype=int)
    jun_end = array(exon_start_list[1:],dtype=int)

    gene_id=ls[0]
    transcript_id = ls[1]

    chr_name = ls[2]

    transcript_candidate_ls = []

    i = 0
    for jun_end in exon_start_list[1:]:
        jun_start = exon_end_list[i]
        jun = chr_name + ":" + jun_start + "_" + jun_end

        if jun_transcript_ID_dict.has_key(jun):
            transcript_candidate_ls.extend(jun_transcript_ID_dict[jun])
        i+=1

    transcript_candidate_set = set(transcript_candidate_ls)



#############################################################################################
    if len(transcript_candidate_set) > 0:
        best_gene_N = 0
        best_gene=""
        for item in transcript_candidate_set:
            if transcript_candidate_ls.count(item)>best_gene_N:
                best_gene_N = transcript_candidate_ls.count(item)
                best_gene = item
        if not name_dt.has_key(transcript_id ):
            name_dt[transcript_id ] = [gene_id,item,["-"]]
    else:
        name_dt[transcript_id ] = [gene_id,"-",["-"]]
target.close()


output = open(output_filename,'w')

for ID in name_dt:
    for refID in name_dt[ID][2]:
        output.write(ID+"\t"+"\t".join(name_dt[ID][0:2])+ "\t"+ refID +"\n")
output.close()

















