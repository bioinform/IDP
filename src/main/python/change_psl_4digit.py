#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 2 :
    psl_filename = sys.argv[1]
else:
    print("usage: psl_filename")
    print("or ")
    sys.exit(1)
################################################################################

psl = open(psl_filename,'r')

for line in psl:
    ls = line.strip().split("\t")
    identity = str(round(float(ls[0])/float(ls[10]),4))
    name = identity + "_" + ls[10]
    if ls[9][-3:] == "ccs":
        name = name + "|ccs" 
    ls[9] = name
    print "\t".join(ls) 
psl.close()
################################################################################
    
