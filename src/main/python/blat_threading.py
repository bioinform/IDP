#!/usr/bin/python
import sys
import struct
import os
import threading
import string

################################################################################

if len(sys.argv) >= 8: #JWDEBUG modified to handle the extra python argument
    blat_threading_pathfilename =  sys.argv[0]
    python_path = sys.argv[1] #JWDEBUG adding this extra agrument to call python
    Nthread1 = int(sys.argv[2])
    blat_option = ' '.join(sys.argv[3:-2])
    SR_pathfilename = sys.argv[-2]
    output_pathfilename = sys.argv[-1]
    
else:
    print("usage: python2.6 blat_threading.py p -t=DNA -q=DNA ~/annotations/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /usr/bin/python intact_SM.fa intact_SM.fa.psl")
    print("or ./blat_threading.py p -t=DNA -q=DNA ~/annotations/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /usr/bin/python intact_SM.fa intact_SM.fa.psl")
    sys.exit(1)
################################################################################
def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    if path == "/":
        path = "./"
    return path, filename
################################################################################
SR_path, SR_filename = GetPathAndName(SR_pathfilename)
output_path, output_filename = GetPathAndName(output_pathfilename)
bin_path2, blat_threading= GetPathAndName(blat_threading_pathfilename)
bin_path1 = bin_path2
################################################################################
SR = open(SR_pathfilename,'r')
SR_NR = 0
for line in SR:
    SR_NR+=1
SR.close()

Nsplitline = 1 + (SR_NR/Nthread1)
if Nsplitline%2==1:
    Nsplitline +=1
ext_ls=[]
j=0
k=0
i=0
while i <Nthread1:
    ext_ls.append( '.' + string.lowercase[j] + string.lowercase[k] )
    k+=1
    if k==26:
        j+=1
        k=0
    i+=1

print "===split SR:==="    
splitSR_cmd = "split -l " + str(Nsplitline) + " " + SR_pathfilename + " " + output_path +SR_filename +"."
print splitSR_cmd
os.system(splitSR_cmd)

##########################################
print "===compress SR.aa:==="    

i=0
T_blat_SR_ls = []
for ext in ext_ls:
    blat_SR_cmd = bin_path1 + "blat " + blat_option + ' ' + output_path + SR_filename + ext + ' ' + output_path + SR_filename + ext + ".psl"
    print blat_SR_cmd
    T_blat_SR_ls.append( threading.Thread(target=os.system, args=(blat_SR_cmd,)) )
    T_blat_SR_ls[i].start()
    i+=1

for T in T_blat_SR_ls:
    T.join()

i=0
T_bestblat_SR_ls = []
for ext in ext_ls:
    bestblat_SR_cmd = python_path + ' ' + bin_path2 + "blat_best.py " + output_path + SR_filename + ext + '.psl 5 > ' + output_path + SR_filename + ext + ".bestpsl"
    print bestblat_SR_cmd
    T_bestblat_SR_ls.append( threading.Thread(target=os.system, args=(bestblat_SR_cmd,)) )
    T_bestblat_SR_ls[i].start()
    i+=1

for T in T_bestblat_SR_ls:
    T.join()

cat_bestpsl_cmd = "cat "
rm_SR_cmd = "rm "
rm_SRpsl_cmd = "rm "
for ext in ext_ls:
    cat_bestpsl_cmd = cat_bestpsl_cmd + output_path + SR_filename + ext + ".bestpsl "
    rm_SR_cmd = rm_SR_cmd + output_path + SR_filename + ext + ' '
    rm_SRpsl_cmd = rm_SRpsl_cmd + output_path + SR_filename + ext + '.psl '
cat_bestpsl_cmd = cat_bestpsl_cmd + " > " + output_pathfilename    
print cat_bestpsl_cmd
print rm_SR_cmd
print rm_SRpsl_cmd
os.system(cat_bestpsl_cmd)
os.system(rm_SR_cmd)
os.system(rm_SRpsl_cmd)
