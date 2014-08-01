#!/usr/bin/python

import sys
import os
import datetime

t0 = datetime.datetime.now()
################################################################################
def compress(seq):
    cps_seq = seq[0]
    pos_ls=[]
    len_ls = []
    r=0
    ref_s = seq[0]
    i=0
    for s in seq:
        if not ref_s == s:
            if r>1:
                len_ls.append(str(r))
                pos_ls.append(str(i)) 
            cps_seq = cps_seq + s
            r=1
            ref_s = s
            i+=1
        else:
            r+=1
    if r>1:
        len_ls.append(str(r))
        pos_ls.append(str(i)) 
    return cps_seq,pos_ls,len_ls            
################################################################################
MinNonN=40
MaxN=1
if len(sys.argv) >= 3:
    for opt in sys.argv[1:]:
        if opt[0]=='-':
            opt_ls = opt.split('=')
            if opt_ls[0]=="-MinNonN":
                MinNonN = int(opt_ls[1])
            elif opt_ls[0]=="-MaxN":
                MaxN = int(opt_ls[1])
    inseq_filename =  sys.argv[-2]
    outseq_prefix = sys.argv[-1]
else:
    print("usage: python compressFASTA.py inseq out_prefix")
    print("or ./compressFASTA.py inseq out_prefix")
    sys.exit(1)

################################################################################

################################################################################
inseq=open(inseq_filename,'r')
outseq=open(outseq_prefix+'cps','w')
idx = open(outseq_prefix+'idx','w')
i=0
I=1
seq=""
for line in inseq:
    if line[0]==">":
        if i>0:
            cps_seq, pos_ls, len_ls = compress(seq)
            NN = cps_seq.count('N')
            if len(cps_seq)-NN>=MinNonN and NN<=MaxN:
                outseq.write(">"+readname+'\n' + cps_seq+'\n')
                idx.write(readname + "\t" + ','.join(pos_ls) + '\t' + ','.join(len_ls) + '\n')
                I+=1
            seq=""
        readname = line[1:-1]
        i+=1
    else:
        seq=seq+line.strip().upper()

print str(datetime.datetime.now()-t0)
cps_seq, pos_ls, len_ls = compress(seq)  
NN = cps_seq.count('N')
if len(cps_seq)-NN>=MinNonN and NN<=MaxN:
    outseq.write(">"+readname+'\n' + cps_seq+'\n')
    idx.write(readname + "\t" + ','.join(pos_ls) + '\t' + ','.join(len_ls) + '\n')

inseq.close()
outseq.close()
idx.close()
print "finsish genome"

################################################################################

