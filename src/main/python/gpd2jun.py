#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 2:
    refFlat_filename = sys.argv[1]
else:
    print("usage: python genephed2gtf.py refFlat.txt source]")
    print("or ./genephed2gtf.py refFlat.txt source")
    sys.exit(1)

################################################################################
    
ref=open(refFlat_filename,'r')
output = open(refFlat_filename+".gtf",'w')
jun_dt = {}
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
    for jun_start in Exon_end_list[:-1]:
        jun_end = Exon_start_list[i]
        output_ls = [chr_name]
        output_ls.append( str(int(jun_start)-50) )
        output_ls.append( str(int(jun_end)+50) )
        output_ls.append(refline_list[0])
        output_ls.append("50")
        output_ls.append(strand)
        output_ls.append( str(int(jun_start)-50) )
        output_ls.append( str(int(jun_end)+50) )
        output_ls.append("0,0,0\t2\t50,50")
        output_ls.append( '0,'+str(int(jun_end)-int(jun_start) +50) )
        locus = chr_name + ':' + jun_start + '_' + jun_end
        jun_dt[locus] = '\t'.join(output_ls)
        i += 1
for locus in jun_dt:
    print jun_dt[locus]
