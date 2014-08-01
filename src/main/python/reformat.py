#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 2:
    filename = sys.argv[1]
else:
    print("usage: ")
    print("or ")
    sys.exit(1)
################################################################################

def replacespacebytab(ls):
    result = []
    for item in ls:
        if not item == "":
            result.append(item)
    return result

################################################################################

file = open(filename,'r')
for line in file:
    if line != "":
        ls = line.strip().split(' ')
        ls = replacespacebytab(ls)
        print '\t'.join(ls)
    else:
        print ""
file.close()
