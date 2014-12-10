#!/usr/bin/python
import sys
import struct
import os
import string
from idpcommon import log_command

################################################################################

if len(sys.argv) >= 9: #JWDEBUG modified to handle the extra python argument
    gmap_threading_pathfilename =  sys.argv[0]
    python_path = sys.argv[1] #JWDEBUG adding this extra agrument to call python
    gmap_path = sys.argv[2]
    gmap_option = ' '.join(sys.argv[3:-3])
    SR_pathfilename = sys.argv[-3]
    gmap_index_pathfoldername = sys.argv[-2]
    output_pathfilename = sys.argv[-1]
    
else:
    print("usage: python2.6 gmap_threading.py -f 1 -D ~/annotations/hg19/UCSC/hg19/Sequence/WholeGenomeFasta/ -d genome/ /usr/bin/python intact_SM.fa intact_SM.fa.psl")
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
bin_path2, gmap_threading= GetPathAndName(gmap_threading_pathfilename)
gmap_folder, gmap_index = GetPathAndName(gmap_index_pathfoldername)
bin_path1 = bin_path2
################################################################################
SR = open(SR_pathfilename,'r')
SR_NR = 0
for line in SR:
    SR_NR+=1
SR.close()

##########################################

gmap_SR_cmd = gmap_path + " " + gmap_option + ' ' + " -D " + gmap_folder + " -d " + gmap_index + " " + SR_pathfilename + ' > ' + output_path + SR_filename + ".psl"
print gmap_SR_cmd
log_command(gmap_SR_cmd)
	
#notice we are not skipping any lines in the PSL file output by gmap (hence the zero in the following command)
bestblat_SR_cmd = python_path + ' ' + bin_path2 + "blat_best.py " + output_path + SR_filename + '.psl 0 > ' + output_pathfilename
log_command(bestblat_SR_cmd)

rm_SRpsl_cmd = "rm " + output_path + SR_filename + '.psl '

print rm_SRpsl_cmd
log_command(rm_SRpsl_cmd)
