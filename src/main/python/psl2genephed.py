#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 3:
    psl_filename = sys.argv[1]
    skip_Nine = int(sys.argv[2])
    Lgap = int(sys.argv[3])
    output_filename = sys.argv[4]

else:
    print("usage:psl2genephed.py psl_file skip_Nine Lgap output_filename")
    print("or ")
    sys.exit(1)
################################################################################
def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    if path == "/":
        path = "./"
    return path, filename

path, filename = GetPathAndName(output_filename)
singleexon_output_filename = path + "singleexon_" + filename
psl = open(psl_filename,'r')
output = open(output_filename,'w')
singleexon_output = open(singleexon_output_filename,'w')
singleexon_psl_output = open(singleexon_output_filename + ".psl",'w')

i=0
for line in psl:
    if i<skip_Nine:
        i+=1
        continue
    result_ls  = []
    jun_start_ls = []
    jun_end_ls = []
    ls = line.strip().split('\t')
    strand =ls[8]
    readname = ls[9]
    leftstart = ls[15]
    rightend = ls[16]
    size_ls = ls[18].strip(',').split(',')
    start_ls = ls[20].strip(',').split(',')

#    print size_ls
#    print start_ls

    chr_name = ls[13]
    j=0
    for start in start_ls[:-1]:

        start = int(start)
        size = int(size_ls[j])
        jun_start = start + size

        jun_end = int(start_ls[j+1])
#        print start
#        print size
#        print jun_start,jun_end
        if jun_end - jun_start >= Lgap:
            jun_start_ls.append(str(jun_start))
            jun_end_ls.append(str(jun_end))
        j+=1

    if 1:
        result_ls.append(readname)
        result_ls.append('.')
        result_ls.append(chr_name)
        result_ls.append(strand)
        result_ls.append(leftstart)
        result_ls.append(rightend)
        result_ls.append(leftstart)
        result_ls.append(rightend)
        result_ls.append(str(len(jun_start_ls)+1))
        jun_end_ls.insert(0,leftstart)
        jun_start_ls.append(rightend)
        result_ls.append(','.join(jun_end_ls)+',')
        result_ls.append(','.join(jun_start_ls)+',')
    if len(jun_start_ls) > 1:
        output.write('\t'.join(result_ls)+'\t'+ ls[0] + '\t' + ls[10] + '\n' )
    else:
        singleexon_output.write('\t'.join(result_ls)+'\t'+ ls[0] + '\t' + ls[10] + '\n' )
        singleexon_psl_output.write(line)
    i+=1
output.close()
psl.close()
singleexon_output.close()
singleexon_psl_output.close()
################################################################################

