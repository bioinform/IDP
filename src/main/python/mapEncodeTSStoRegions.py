#!/usr/bin/python
import sys
from operator import itemgetter, attrgetter
import bisect
import re

TSS_THRESHOLD = 0.99
OVERLAP_THRESHOLD = 0.0

### parse_tss_file
##################
def parse_tss_file(tss_file):
    
    tss_dict = {}
        
    for line in tss_file:
        fields = line.split()
        rname = fields[0]

        start_pos = int(fields[1])
        end_pos = int(fields[2]) - 1   # included
        #description = fields[3]
        #exp_level = fields[6]
        #idr = fields[7]
        if not tss_dict.has_key(rname):
            tss_dict[rname] = []
        tss_dict[rname].append([start_pos, end_pos, line[:-1]])
        
    return [tss_dict] 


### sort_tss_dict
##################
def sort_tss_dict(tss_dict):
    
    tss_start_pos_dict = {}
    tss_start_sorted_index_dict = {}
    tss_end_pos_dict = {}
    tss_end_sorted_index_dict = {}
    
    CHR_LIST = tss_dict.keys()
     
    for chr in CHR_LIST:     
        # Sort based on start position
        temp_list = [tss_dict[chr][j] + [j] for j in range(len(tss_dict[chr]))]
        sorted_temp_list = sorted(temp_list, key=itemgetter(0))
        tss_start_pos_dict[chr] = [sorted_temp_list[j][0] for j in range(len(sorted_temp_list))]
        tss_start_sorted_index_dict[chr] = [sorted_temp_list[j][-1] for j in range(len(sorted_temp_list))]
        # Sort based on end position
        sorted_temp_list = sorted(temp_list, key=itemgetter(1))
        tss_end_pos_dict[chr] = [sorted_temp_list[j][1] for j in range(len(sorted_temp_list))]
        tss_end_sorted_index_dict[chr] = [sorted_temp_list[j][-1] for j in range(len(sorted_temp_list))]

    return [tss_start_pos_dict, tss_start_sorted_index_dict, tss_end_pos_dict, tss_end_sorted_index_dict] 



### parse_region_file
#####################
def parse_region_file(region_file):
    
    region_dict = {}
        
    for line in region_file:
        fields = line.split()
        rname = re.split(r'_|:', fields[0])[0]

        description = fields[0]
        regions = [int(i) for i in fields[1:]]
        
        regions_list = []
        for i in range(len(regions) - 1):
            regions_list.append([regions[i], regions[i+1] - 1]) 
        
        if not region_dict.has_key(rname):
            region_dict[rname] = []
        region_dict[rname].append([regions_list] + [description])
        
    return [region_dict]


### overlap
############
def overlap(region, tss):
    
    start_point = max(region[0], tss[0])
    end_point = min(region[1], tss[1])
    
    if (end_point < start_point):
         return 0
    
    return (1. * (end_point - start_point) / (tss[1] - tss[0])) 
    


### map_tss_to_regions
########################
def map_tss_to_regions(tss_dict, region_dict, tss_start_pos_dict, tss_start_sorted_index_dict, 
                                              tss_end_pos_dict, tss_end_sorted_index_dict, overlap_threshold ): 
    output_dict = {}
    tss_dict_len = len(tss_start_pos_dict);
    
    CHR_LIST = region_dict.keys()
    for chr in CHR_LIST:
        
        output_dict[chr] = {}
        for gene in region_dict[chr]:
            #print '***********'
            output_dict[chr][gene[1]] = {} 
            
            if (tss_start_pos_dict.has_key(chr)):
                start_index = bisect.bisect_right(tss_start_pos_dict[chr], gene[0][-1][1])
                end_index = bisect.bisect_left(tss_end_pos_dict[chr], gene[0][0][0])
                tss_candidates = set(tss_start_sorted_index_dict[chr][:start_index]) & set(tss_end_sorted_index_dict[chr][end_index:])
            else:
                tss_candidates = set()
            #print tss_candidates
            #print gene[0]

            for region in gene[0]: 
                region_name = str(region[0]) + '-' + str(region[1])
                output_dict[chr][gene[1]][region_name] = []
                for tss_idx in tss_candidates:
                    if (overlap(region, tss_dict[chr][tss_idx][0:2]) >  overlap_threshold):
                        tss = tss_dict[chr][tss_idx][-1] + '\t' + str(overlap(region, tss_dict[chr][tss_idx][0:2])) + '\n'
                        output_dict[chr][gene[1]][region_name].append(tss)
                        #print gene[1]
                        #print region, overlap(region, tss_dict[chr][tss_idx][0:2])
                        #print tss_dict[chr][tss_idx][0:2]
                        #print '-----------------------'
                    
        
    return output_dict


###generate_output_file
#########################
def filter_mapped_tss_prediction(output_dict, threshold):
    
    CHR_LIST = output_dict.keys()
    for chr in CHR_LIST:
        for gname in sorted(output_dict[chr].keys()):
            gene = output_dict[chr][gname]
            for region_name in gene.keys():
                if (len(gene[region_name]) == 0):
                    continue

                region_tss = []
                for tss in gene[region_name]:
                    fields = tss.split()
                    idr = float(fields[7])
                    tss_pred_strength = float(fields[3].split(':')[5])
                    if (tss_pred_strength >=  threshold):
                        region_tss.append(tss)
                        
                gene[region_name] = region_tss

    return output_dict


###generate_output_file
#########################
def generate_output_file(output_file_str, output_dict):
    
    output_file = open(output_file_str, 'w')
    CHR_LIST = output_dict.keys()
    for chr in CHR_LIST:
        for gname in sorted(output_dict[chr].keys()):
            gene = output_dict[chr][gname]
            gene_keys = sorted(gene.keys())
            
            output_file.write(gname + '\t' + chr + '\t' + str(len(gene.keys())) + '\n')
            
            for region_name in gene_keys:
                output_file.write(region_name.ljust(30))
            output_file.write('\n')
            for region_name in gene_keys:
                output_file.write(str(len(gene[region_name])).ljust(30))
            output_file.write('\n')
            for region_name in gene_keys:
                for tss in gene[region_name]:
                    output_file.write(tss)
    
    output_file.close()

### Main
########
def main():
    # Read input parameters
    tss_file_str = sys.argv[1]
    region_file_str = sys.argv[2]
    
    output_file_str = 'encodeTSS_mapped_regions.txt'
    
    tss_file = open(tss_file_str, 'r')
    [tss_dict] = parse_tss_file(tss_file)
    [tss_start_pos_dict, tss_start_sorted_index_dict, tss_end_pos_dict, tss_end_sorted_index_dict]  = sort_tss_dict(tss_dict)
    tss_file.close()
    
    region_file = open(region_file_str, 'r')
    [region_dict] = parse_region_file(region_file)
    region_file.close()
    
    output_dict = map_tss_to_regions(tss_dict, region_dict, tss_start_pos_dict, tss_start_sorted_index_dict, 
                                     tss_end_pos_dict, tss_end_sorted_index_dict, OVERLAP_THRESHOLD)
    output_dict = filter_mapped_tss_prediction(output_dict, TSS_THRESHOLD)
    generate_output_file(output_file_str, output_dict)


if __name__ == '__main__':
    main()
