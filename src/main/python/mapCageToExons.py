#!/usr/bin/python
import sys
from operator import itemgetter, attrgetter
import bisect
import re
import os
from binaidp import log_command

CAGE_margin = 100


### Main
########
def main():
    # Read input parameters
    tss_file_str = sys.argv[1]
    exon_file_str = sys.argv[2]
    ref_gpd_str = sys.argv[3]
    I_ref5end_isoformconstruction = int(sys.argv[4])
    output_file_str = sys.argv[5]

    
    exon_file = open(exon_file_str, 'r')
    exon_bed_file = open(exon_file_str + ".map2cage.bed", 'w')
    for line in exon_file:
        fields = line.strip().split("\t")
        chrname = fields[0].split(":")[0]
        for idx in range(1, len(fields)):
            exon_bed_file.write('\t'.join([chrname, str(int(fields[idx]) - CAGE_margin), str(int(fields[idx]) + CAGE_margin), fields[0]]) + '\n')
    exon_bed_file.close()
    exon_file.close()

    bedtools_cmnd = "bedtools intersect -wb -a " + tss_file_str + " -b " + exon_file_str + ".map2cage.bed > " + output_file_str + ".all"
    os.system(bedtools_cmnd)

    # Remove redundancy
    if (I_ref5end_isoformconstruction):
        ref_file = open(ref_gpd_str, 'r')
        ref_bed_file = open(ref_gpd_str + ".map2cage.bed", 'w')
        for line in ref_file:
            fields = line.strip().split("\t")
            if (fields[3] == "+"):
                idx = 4
            else:
                idx = 5
            chrname = fields[2]
            ref_bed_file.write('\t'.join([chrname, str(int(fields[idx]) - CAGE_margin), str(int(fields[idx]) + CAGE_margin), fields[0]]) + '\n')
        ref_bed_file.close()
        ref_file.close()
        
        bedtools_cmnd = "bedtools intersect -v -a " + output_file_str + ".all" + " -b " + ref_gpd_str + ".map2cage.bed > " + output_file_str 

        os.system(bedtools_cmnd)


if __name__ == '__main__':
    main()
