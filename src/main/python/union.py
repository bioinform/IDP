#!/usr/bin/python

import sys
import os
import datetime
import math

#############################################################################################
def Genelist2Set(filename):
    result=set()
    file=open(filename,'r')
    temp_set=set(file)
    for item in temp_set:
        #result.add(item.strip().strip('"').upper())
        result.add(item.strip().strip('"'))
    file.close()
    return result

if len(sys.argv)>=3:
    filename1=sys.argv[1]
    filename2=sys.argv[2]

else:
    print("usage: ./intersection.py filename1 filename2")
    sys.exit(1)

set1=Genelist2Set(filename1)
set2=Genelist2Set(filename2)
result_set= set1 | set2
for item in result_set:
    print item +'\t'+ filename1
