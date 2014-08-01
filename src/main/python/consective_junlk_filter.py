#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 4 :
    ref_gpd_filename = sys.argv[1]
    tag_gpd_filename = sys.argv[2]
    output_filename_suffix = sys.argv[3]
else:
    print("usage:ref_gpd_filename tag_gpd_filename output_filename_suffix")
    print("or ./")
    sys.exit(1)

def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])
    if path == "":
        path = "."
    path = path +'/'
    return path, filename

output_filename_suffix_path, output_filename_suffix_filename = GetPathAndName(output_filename_suffix)
novel_output_filename = output_filename_suffix_path + "negative_" + output_filename_suffix_filename 
known_output_filename = output_filename_suffix_path + "maypositive_" + output_filename_suffix_filename

################################################################################
print "loading reference annotation", ref_gpd_filename    
ref=open(ref_gpd_filename,'r')
ref_junlk_dt={}
for refline in ref:
    refline_list=refline.split()
    exon_start_list=refline_list[9].strip(',').split(',')
    exon_start_list=exon_start_list[1:]
    exon_end_list=refline_list[10].strip(',').split(',')
    exon_end_list=exon_end_list[:-1]

    if not ref_junlk_dt.has_key(refline_list[2]):
        ref_junlk_dt[refline_list[2]]={}

    i=1
    while i < len(exon_end_list):
        start = exon_end_list[i]
        end=exon_start_list[i]
        next_start = exon_end_list[i]
        next_end = exon_start_list[i]
        junlk_str = str(start) + "_" +  str(end) + "_" + str(next_start) + "_" + str(next_end) 
        if not ref_junlk_dt[refline_list[2]].has_key(junlk_str):
            ref_junlk_dt[refline_list[2]][junlk_str] = 0
        ref_junlk_dt[refline_list[2]][junlk_str] += 1
        i+=1
            
ref.close()

################################################################################
def check_novelity(ref_junlk_dt,chr_name,exon_start_list,exon_end_list):
#    0 is novel, at least one of the twojun_lk is novel 
#    1 is known, all twojun_lk is known
    i=1
    while i < len(exon_end_list):
        start = exon_end_list[i]
        end=exon_start_list[i]
        next_start = exon_end_list[i]
        next_end = exon_start_list[i]
        junlk_str = str(start) + "_" +  str(end) + "_" + str(next_start) + "_" + str(next_end)
        if not ref_junlk_dt[chr_name].has_key(junlk_str):
            return 0
        i+=1
    return 1
################################################################################
novel_output = open(novel_output_filename,'w')
known_output = open(known_output_filename,'w')
print "loading target annotation", tag_gpd_filename

tag = open(tag_gpd_filename,'r')
for line in tag:
    line_list=line.split()
    exon_start_list=line_list[9].strip(',').split(',')
    exon_start_list=exon_start_list[1:]
    exon_end_list=line_list[10].strip(',').split(',')
    exon_end_list=exon_end_list[:-1]

    if len(exon_end_list) == 1:
        print "single-junction transcript:", line.strip()
        continue
    elif len(exon_end_list) == 0:
        print "single-exon transcript:", line.strip()
        continue
    if not ref_junlk_dt.has_key(line_list[2]):
        print "no reference chromosome", line_list[2]
        continue #novel

    I = check_novelity(ref_junlk_dt,line_list[2],exon_start_list,exon_end_list)

    if I == 0:
        novel_output.write(line)
    elif I == 1:
        known_output.write(line)

tag.close()
novel_output.close()
known_output.close()
