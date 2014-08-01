#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 3:
    ref_refFlat_filename = sys.argv[1]
    tag_refFlat_filename = sys.argv[2]
    output_filename = sys.argv[3]
else:
    print ("usage: python novel_genephed.py ref_refFlat_filename tag_refFlat_filename  output_filename")
    print ("or ./novel_genephed.py ref_refFlat_filename tag_refFlat_filename output_filename")
    sys.exit(1)

################################################################################
def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    if len(ls)==1:
        path="./"
    return path, filename

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
addrefFlat(ref_iso_dt,ref_dt, ref_refFlat_filename)
################################################################################

output = open(output_filename,'w')
tag = open(tag_refFlat_filename,'r')
known_locus_set = set()
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

    old_end = ls[10].strip(',').split(',')[-1]
    old_start = ls[9].strip(',').split(',')[0]


    if ref_dt.has_key(locus):
        known_locus_set.add(locus)
    if ref_iso_dt.has_key(locus):
        locus_ls = locus.split(':')
        output_ls = []
        output_ls.append('|'.join(ref_iso_dt[locus][0]))
        output_ls.append('|'.join(ref_iso_dt[locus][1]))
        output_ls.append(locus_ls[0])
        output_ls.append('?')
        ls1ls =locus_ls[1].split("-")
        exon_start_ls=ls1ls[1].split('_')
        exon_end_ls=ls1ls[0].split('_')




        oldexon_start_ls = exon_start_ls
        oldexon_end_ls = exon_end_ls
        oldexon_start_ls.insert(0,old_start)
        oldexon_end_ls.append(old_end )

        temp_output_ls = output_ls
        N_exon = str( len(exon_end_ls) )
        temp_output_ls.append( exon_start_ls[0] )
        temp_output_ls.append( exon_end_ls[-1] )
        temp_output_ls.append( exon_start_ls[0] )
        temp_output_ls.append( exon_end_ls[-1] )
        temp_output_ls.append( N_exon )

        temp_output_ls.append( ','.join(oldexon_start_ls) + ',' )
        temp_output_ls.append( ','.join(oldexon_end_ls) + ',' )
        temp_output_ls[3] = "+"
        print '\t'.join(temp_output_ls) + '\t' + gene + '\t'.join(ls[10:])



        exon_start_ls.insert(0, str( int(exon_end_ls[0])-100 ) )
        exon_end_ls.append( str( int(exon_start_ls[-1])+100 ) )

        N_exon = str( len(exon_end_ls) )
        output_ls.append( exon_start_ls[0] )
        output_ls.append( exon_end_ls[-1] )
        output_ls.append( exon_start_ls[0] )
        output_ls.append( exon_end_ls[-1] )
        output_ls.append( N_exon )
        output_ls.append( ','.join(exon_start_ls) + ',')
        output_ls.append( ','.join(exon_end_ls) + ',')

#       print '\t'.join(temp_output_ls) + '\t' + gene + '\t' + ls[11] + '\t' + ls[12]

        continue
    output.write( line  )
tag.close()
output.close()

ref_refFlat_path, ref_refFlat_file = GetPathAndName(ref_refFlat_filename)

tag_refFlat_path, tag_refFlat_file = GetPathAndName(tag_refFlat_filename)
known_output = open(tag_refFlat_path + "known_"+ tag_refFlat_file + '_' + ref_refFlat_file,'w')
for locus in known_locus_set:
    for refline in ref_dt[locus]:
        known_output.write(refline)
known_output.close()
