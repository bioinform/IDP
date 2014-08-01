#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 3:
    filename = sys.argv[1]
    pythonpath = sys.argv[2]
else:
    print("usage: ./changetitle.py filename pythonpath")
    print("or python changetitle.py filename pythonpath")
    sys.exit(1)
################################################################################
ls = []
file = open(filename,'r')
i=0
for line in file:
    if i == 0:
        ls.append( pythonpath + "\n" )
        i+=1
        continue
    ls.append( line )

file.close()
################################################################################
file = open(filename,'w')
for item in ls:
    file.write(item)
file.close()
