#!/usr/bin/python

from numpy import *
import sys
import os
from copy import *

if len(sys.argv) >= 3:
    output_filename = sys.argv[1]
    iso_candidate_filename_ls = sys.argv[2:]
else:
    print("usage: python merge_isoform_construction.py [iso_candidate_iso.aa]")
    print("or ./merge_isoform_construction.py [iso_candidate_iso.aa]")
    sys.exit(1)

################################################################################################
candidate_dt={}
for iso_candidate_filename in iso_candidate_filename_ls: 
    iso_candidate_file = open(iso_candidate_filename,'r')
    for line in iso_candidate_file:
        ls = line.strip().split("kinfai")
        locus = ls[0]
        item = ls[1]
        gpd_strand = ls[2]
        twoending_str = ls[3]

        if not candidate_dt.has_key(locus):
            candidate_dt[locus] = {}
        if not candidate_dt[locus].has_key(item):
            candidate_dt[locus][item] ={}
        if not candidate_dt[locus][item].has_key(gpd_strand):
            candidate_dt[locus][item][gpd_strand] = set()
        candidate_dt[locus][item][gpd_strand].add(twoending_str)

    iso_candidate_file.close()

################################################################################################
testout = open(output_filename,'w')
for locus in candidate_dt:

#    for item in candidate_dt[locus]:
#        for gpd_strand in candidate_dt[locus][item]:
#            for twoending_str in candidate_dt[locus][item][gpd_strand]:
#                print gpd_strand + '\t' + twoending_str +'\t'+ item



    i = 1
    locus_ls = locus.split(":")
    chr_name = locus_ls[0]
    for item in candidate_dt[locus]:
        temp = [locus]
        temp.append(chr_name)

        ls  = item.strip("_").split('_')
        L_ls = len(ls)/2
        exon_start_ls = []
        exon_end_ls = []

        k = 0
        while k < L_ls:
            exon_end_ls.append(  ls[2*k] )
            exon_start_ls.append( ls[2*k+1] )
            k+=1


        for gpd_strand in candidate_dt[locus][item]:
            if gpd_strand == '':
                temp.append("+")
            else:
                temp.append(gpd_strand)               
            for twoending_str in candidate_dt[locus][item][gpd_strand]:
                first_temp =deepcopy(temp)
                first_temp.insert(1,  locus + '.' +str(i) )

                twoending_ls = twoending_str.split('_')
                if twoending_ls[0] == "-1" and twoending_ls[1] == "-1":
                    fristexonstart = str(  int(exon_end_ls[0])-200  )
                    lastexonend =  str(  int(exon_start_ls[-1])+200  )
                elif twoending_ls[0] != "-1" and twoending_ls[1] == "-1":
                    fristexonstart = twoending_ls[0] 
                    lastexonend =  str(  int(exon_start_ls[-1])+200  )                 
                elif twoending_ls[0] == "-1" and twoending_ls[1] != "-1":
                    fristexonstart = str(  int(exon_end_ls[0])-200  )
                    lastexonend =  twoending_ls[1]                
                elif twoending_ls[0] != "-1" and twoending_ls[1] != "-1":
#                    print "two ending: ", item
                    fristexonstart = twoending_ls[0] 
                    lastexonend =  twoending_ls[1]                
                temp_exon_start_ls = deepcopy (exon_start_ls )
                temp_exon_end_ls = deepcopy(exon_end_ls)
                temp_exon_start_ls.insert(0,fristexonstart ) 
                temp_exon_end_ls.append( lastexonend )
                second_temp = []
                second_temp.append(fristexonstart)
                second_temp.append(lastexonend)
                second_temp.append(fristexonstart)
                second_temp.append(lastexonend)
                second_temp.append( str(L_ls+1) )
                second_temp.append( ','.join(temp_exon_start_ls) + ',' )
                second_temp.append( ','.join(temp_exon_end_ls) + ',' )

                testout.write( '\t'.join(first_temp) + '\t' + '\t'.join(second_temp) + '\n' )
                i +=1
       
testout.close()        
#    exon_boundary_ls = exon_dt[locus]
#    LR_startend_ls = searchsorted(exon_boundary_ls,[LR_start,LR_end])
#    LR_start_index = LR_startend_ls[0]
#    LR_end_index = LR_startend_ls[1]
#    start_end_lk_ls,start_end_lk_ls = construct_jun_linkage(jun_start_dt[chr_name],jun_end_dt[chr_name],gene_start,gene_end,exon_boundary_ls,LR_start,LR_end,LR_start_index,LR_end_index,start_terminus_s,end_terminus_s)


################################################################################

