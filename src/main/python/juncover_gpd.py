#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 4 :
    ref_jun_filename = sys.argv[1]
    tag_gpd_filename = sys.argv[2]
    output_filename = sys.argv[3]
else:
    print("usage: python ~/3seq/juncover.py ref_jun.bed target.gpd")
    print("or ~/3seq/juncover.py ref_jun.bed target.gpd")
    sys.exit(1)

################################################################################

ref_jun=open(ref_jun_filename,'r')
ref_jun_dt = {}

for line in ref_jun:
    if line[0:5]=='track':
        continue
    else:
        line_list=line.strip().split("\t")
        chr_name = line_list[0]
        thickness=line_list[10].split(',')
        leftpos=str(int(line_list[1])+int(thickness[0]))
        rightpos=str(int(line_list[2])-int(thickness[1]))
        locus = chr_name + ":" + leftpos + "_" +rightpos
        if not ref_jun_dt.has_key(locus):
            ref_jun_dt[locus] = 0
        ref_jun_dt[locus] += 1
ref_jun.close()

################################################################################

gpd = open(tag_gpd_filename,'r')
output = open(output_filename + ".bed",'w')
for refline in gpd:
    refline_list=refline.split()
    Exon_start=int(refline_list[4])
    Exon_end=int(refline_list[5])
    Exon_start_list=refline_list[9].strip(",").split(',')
    Exon_end_list=refline_list[10].strip(",").split(',')
    strand=refline_list[3]
    chr_name = refline_list[2]
    j = 1 
    jun_start_ls = []
    jun_end_ls = []
    for jun_start in Exon_end_list[:-1]:
        jun_start_ls.append(jun_start)
        jun_end = Exon_start_list[j]
        jun_end_ls.append(jun_end)
        j += 1

    j=0
    I = 0
    for start in jun_start_ls:
        end=jun_end_ls[j]
        locus = chr_name + ":" + start + "_" + end
        if ref_jun_dt.has_key(locus):
            I += 1
        else:
            jun_start = int(start)
            jun_end = int(end)
            output_ls = [chr_name,str(jun_start-50),str(jun_end+50),"LR_covered_jun","50","+",str(jun_start-50),str(jun_end+50),"0,0,0","2","50,50"]
            jun_L = jun_end - jun_start - 50
            output_ls.append("0,"+str(jun_L))
            output.write( "\t".join(output_ls) +'\n')
              
        j+=1
    if I == len(jun_start_ls):
        print refline.strip()
    
gpd.close()
################################################################################
output.close()
