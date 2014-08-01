#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 3 :
    psl_filename = sys.argv[1]
    threeend_filename =sys.argv[2]

else:
    print("usage: psl_filename threeend_filename")
    print("or ")
    sys.exit(1)
################################################################################
threeend_dt = {}
threeend = open(threeend_filename,'r')
for line in threeend:
    ls = line.strip().split('\t')
    readname = ls[0]
    strand = ls[1]
    if threeend_dt.has_key(readname):
        print "error: two 3end of\t" + readname
    threeend_dt[readname] = strand
threeend.close()
################################################################################
psl = open(psl_filename,'r')

for line in psl:
    ls = line.strip().split("\t")
    identity = str(round(float(ls[0])/float(ls[10]),4))
    name = identity + "_" + ls[10]
    strand = ls[8]
    readname = ls[9]
    size_ls = ls[18].strip(',').split(',')
    readstart_ls = ls[19].strip(',').split(',')
    start_ls = ls[20].strip(',').split(',')
    end_pt = int(readstart_ls[-1]) + int(size_ls[-1]) 

    right_L = int(ls[10]) - int(ls[12])
    left_L = int(ls[11])


    last_ls = ls[9].split("+")
    if len(last_ls) == 2:
        if last_ls[0] == "ccs":
            Iccs = 1
        else:
            Icss = 0
    else:
        if ls[9][-3:] == "ccs":
            Iccs = 1
        else:
            Icss = 0


    if ls[9][-3:] == "ccs":
        name = name + "|ccs"

    if threeend_dt.has_key(readname):
#        print name, threeend_dt[readname], strand
        if threeend_dt[readname] == strand:
            name = name+"+"+str(right_L)
        else:
            name = name+"-"+str(left_L)

    ls[9] = name
    print "\t".join(ls) 
psl.close()
################################################################################
    
