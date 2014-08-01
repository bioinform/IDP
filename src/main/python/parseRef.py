#!/usr/bin/python
from operator import itemgetter, attrgetter
from datetime import *
import sys


"""
Note: 
"""

### Compute number of bases to junction gap from an exon index (exclusive)
### exon_list format:  [start+point, end_point, len]
##########
def comp_len_to_junc_gap_forward(start_index, exon_list, isoform_exons, isoform_points_dict):
    if (start_index >= len(exon_list)):
        print 'Warning: Called comp_len_to_junc_gap_forward with start_idx of out of bound' 
        return 0
    index = start_index
    l = 0
    next_index = index +1
    while (next_index< len(exon_list)):
        region_name_temp = 'P' + str(isoform_points_dict[exon_list[next_index][0]]) + ':P' + str(isoform_points_dict[exon_list[next_index][1]])
        if (region_name_temp in isoform_exons):   # this exon is in the isoform
            if ((exon_list[next_index][0] - 1) == exon_list[index][0]):
                l += exon_list[next_index][2]
                index = next_index
            else:
                break
        next_index += 1
    return l
##########
def comp_len_to_junc_gap_backward(start_index, exon_list, isoform_exons, isoform_points_dict):
    if (start_index >= len(exon_list)):
        print 'Warning: Called comp_len_to_junc_gap_backward with start_idx of out of bound' 
        return 0
    index = start_index
    l = 0
    next_index = index -1
    while (next_index >= 0):
        region_name_temp = 'P' + str(isoform_points_dict[exon_list[next_index][0]]) + ':P' + str(isoform_points_dict[exon_list[next_index][1]])
        if (region_name_temp in isoform_exons):   # this exon is in the isoform
            if ((exon_list[next_index][1] + 1) == exon_list[index][0]):
                l += exon_list[next_index][2]
                index = next_index
            else:
                break
        next_index -= 1
    return l
##########

### 
##########
def compute_min_ljust(str_value):
    return (len(str_value)/10 + 1) * 10


### Generate gene exons
### Note: RefSeq has range (start_point, end_point]
### Note: RefSeq exon file  has range [start_point, end_point]
##########
def generate_gene_exons(ref_file_str, exon_file_str):

    file_read = open(ref_file_str, 'r')
    file_write = open(exon_file_str, 'w')
 
    # Keep exon regions for each genome
    gene_dict = dict()

    for line in file_read:
        fields = line.split()
        rname = fields[2]
        gname = fields[0]
        num_exons = int(fields[8])
        start_pos = fields[9].split(',')[:-1]  # The sequence ends in ,
        end_pos = fields[10].split(',')[:-1]

        exon_list = []
        for i in range(num_exons):
            exon_list.append([int(start_pos[i]) + 1, int(end_pos[i])])
        if (not gene_dict.has_key(rname)):
            gene_dict[rname] = dict()
        if gene_dict[rname].has_key(gname):
            exon_list += gene_dict[rname][gname]
            # print gname + ' repeated.'
        gene_dict[rname][gname] = exon_list  # Save the list of exons for each gene

    for chr in gene_dict.keys():
        for gname in gene_dict[chr].keys():
            exon_sorted = sorted(gene_dict[chr][gname], key=itemgetter(0, 1))  # Sort by start position then end position
            i = 1
            while (i < len(exon_sorted)):
                if ((exon_sorted[i][0] == exon_sorted[i - 1][0]) and 
                    (exon_sorted[i][1] == exon_sorted[i - 1][1])):
                    del exon_sorted[i]  # Delete the repeated exon
                elif (exon_sorted[i][0] <= exon_sorted[i - 1][1]):
                    temp_exons = sorted([exon_sorted[i][0], exon_sorted[i - 1][0], exon_sorted[i][1], exon_sorted[i - 1][1]])
                    del exon_sorted[i - 1]  # Delete the two overlapping exons
                    del exon_sorted[i - 1]
                    if (temp_exons[0] == temp_exons[1]):  # Based on two exons overlap type, re-generate 2 or 3 new ones
                        exon_sorted.insert(i - 1, [temp_exons[1], temp_exons[2]])
                        exon_sorted.insert(i , [temp_exons[2] + 1, temp_exons[3]])
                    elif (temp_exons[2] == temp_exons[3]):
                        exon_sorted.insert(i - 1, [temp_exons[0], temp_exons[1] - 1])
                        exon_sorted.insert(i , [temp_exons[1], temp_exons[2]])
                    else:
                        exon_sorted.insert(i - 1, [temp_exons[0], temp_exons[1] - 1])
                        exon_sorted.insert(i, [temp_exons[1], temp_exons[2]])
                        exon_sorted.insert(i + 1, [temp_exons[2] + 1, temp_exons[3]])
                    exon_sorted = sorted(exon_sorted, key=itemgetter(0, 1))  # re-sort the exon positions
                else:
                    i += 1  # After sorting re-evaluate the same index unless there is no change
            gene_dict[chr][gname] = exon_sorted

    keys_sorted = sorted(gene_dict.keys())
    
    # Print out the exon positions sorted in chr name and then gene start position
    for chr in keys_sorted:
        gene_list = []
        for gname in gene_dict[chr].keys():
            exon_sorted = sorted(gene_dict[chr][gname], key=itemgetter(0, 1))  # Sort by start position (There should not be any overlap between exons)
            gene_dict[chr][gname] = exon_sorted
            gene_list.append([gname, gene_dict[chr][gname][0][0]])
        keys = sorted(gene_list, key=itemgetter(1)) 
        for gname in keys:
            exon_sorted = sorted(gene_dict[chr][gname[0]], key=itemgetter(0, 1))  # Sort by start position (There should not be any overlap between exons)
            gene_dict[chr][gname[0]] = exon_sorted
            if (len(gname[0]) < 20):
                file_write.write(gname[0].ljust(20))
            else:
                file_write.write(gname[0].ljust((len(gname[0]) / 10 + 1) * 10))
            file_write.write(str(len(exon_sorted)).ljust(20) + chr.ljust(20) + '\n')
            for i in range(len(exon_sorted)):
                file_write.write(str(exon_sorted[i][0]).ljust(20) + str(exon_sorted[i][1]).ljust(20) + (str(exon_sorted[i][1] - exon_sorted[i][0] + 1).ljust(20)))
                file_write.write('\n')

    file_read.close()
    file_write.close() 

# #########


###
###################
def sanity_check_isoform_regions_length(gene_exons_dict, gene_regions_dict, gene_isoforms_dict, 
                                        isoforms_regions_len_dict, genes_regions_len_dict):
    
    # Sanity check the isoforms regions length
    for i in gene_exons_dict.keys():
        for j in gene_exons_dict[i].keys():
            for k in gene_regions_dict[i][j]:
                if not isoforms_regions_len_dict[i][j].has_key(k):
                    continue
                region_len = 0
                for m in gene_isoforms_dict[i][j]:
                    if not isoforms_regions_len_dict[i][j][k].has_key(m):
                        continue
                    if (region_len == 0):
                        if (isoforms_regions_len_dict[i][j][k][m] > 0):
                            genes_regions_len_dict[i][j][k] = isoforms_regions_len_dict[i][j][k][m]
                            region_len = genes_regions_len_dict[i][j][k]
                    elif (isoforms_regions_len_dict[i][j][k][m] > 0):
                        if (genes_regions_len_dict[i][j][k] != isoforms_regions_len_dict[i][j][k][m]):
                            print 'Isoforms do not match in region length, gene:  %s, chr %s' % (j, i)
                            exit(1)
#######################

###
#############
def compute_isoform_length(start_pos, end_pos):
    length = 0;

    for i in range(len(start_pos)):
        length += (end_pos[i] - start_pos[i])   # In refSeq format, end point is not included
        
    return length
##########


# ## Generate gene isoforms
# #########
def generate_gene_regions(ref_file_str, exon_file_str, READ_LEN, READ_JUNC_MIN_MAP_LEN):

    file_read = open(exon_file_str, 'r')

    # Keep the list of exons for each gene 
    gene_exons_dict = dict()
    gene_points_dict = dict()
    gene_regions_dict = dict()
    isoforms_regions_len_dict = dict()
    genes_regions_len_dict = dict()
        
    num_exons = 0
    for line in file_read: 
        if (num_exons == 0): 
            fields = line.split()
            gname = fields[0]
            rname = fields[2]
            num_exons = int(fields[1])
            exon_list = []
            points_dict = dict()
            point_index = 0
            if (not gene_regions_dict.has_key(rname)):
                gene_regions_dict[rname] = dict()
                isoforms_regions_len_dict[rname] = dict()
                genes_regions_len_dict[rname] = dict()
                
            gene_regions_dict[rname][gname] = dict()
            isoforms_regions_len_dict[rname][gname] = dict()
            genes_regions_len_dict[rname][gname] = dict()
        else:
            fields = line.split()
            exon_list.append([int(fields[0]), int(fields[1]), int(fields[2])])
            points_dict[int(fields[0])] = point_index
            point_index += 1
            points_dict[int(fields[1])] = point_index
            point_index += 1
            if (not gene_exons_dict.has_key(rname)):
                gene_exons_dict[rname] = dict()
                gene_points_dict[rname] = dict()
            gene_exons_dict[rname][gname] = exon_list
            gene_points_dict[rname][gname] = points_dict
            num_exons -= 1;

    file_read.close()

    # Define input/output file
    file_read = open(ref_file_str, 'r')

    # Keep isoforms per gene and exons per isoform 
    gene_isoforms_dict = dict()
    gene_isoforms_length_dict = dict()
        
    for line in file_read:
        fields = line.split()
        rname = fields[2]
        gname = fields[0]
        isoform_name = fields[1]
        num_exons = int(fields[8])
        start_pos = [int(x) for x in fields[9].split(',')[:-1]]  # The sequence ends in ,
        end_pos = [int(x) for x in fields[10].split(',')[:-1]]
        # Generate list of isoforms
        isoform_list = [isoform_name]
        if (not gene_isoforms_dict.has_key(rname)):
            gene_isoforms_dict[rname] = dict()
            gene_isoforms_length_dict[rname] = dict()
        if gene_isoforms_dict[rname].has_key(gname):
            isoform_list += gene_isoforms_dict[rname][gname]
        gene_isoforms_dict[rname][gname] = isoform_list
        gene_isoforms_length_dict[rname][gname + '_' + isoform_name] = compute_isoform_length(start_pos, end_pos) 
        
        # Generate exon indicator for the isoform
        exon_list = gene_exons_dict[rname][gname]  # Note: Assuming sorted start/end positions 
        region_ind = gene_regions_dict[rname][gname]
        isoform_exons = []
        
        # Check the exon regions
        j = 0
        for i in range(len(exon_list)):
            flag = (j < num_exons)
            while (flag):
                p0 = exon_list[i][0]
                p1 = exon_list[i][1]
                l = exon_list[i][2]
                if ((p0 >= start_pos[j]) and
                    (p1 <= end_pos[j])):
                    region_name = 'P' + str(gene_points_dict[rname][gname][p0]) + ':' + 'P' + str(gene_points_dict[rname][gname][p1])
                    if (l >= READ_LEN):  # A valid region for exon
                        temp_isoform_name = set()
                        temp_isoform_name.add(isoform_name)
                        #temp_isoform_name = {isoform_name}
                        if (region_ind.has_key(region_name)):
                            temp_isoform_name = temp_isoform_name.union(region_ind[region_name])
                        else:
                            isoforms_regions_len_dict[rname][gname][region_name] = dict()
                        isoforms_regions_len_dict[rname][gname][region_name][isoform_name] = l - READ_LEN + 1
                        region_ind[region_name] = temp_isoform_name
                    isoform_exons.append(region_name)
                    flag = False
                elif (p0 == start_pos[j]):
                    print 'Out-of-order exon position for isoform ' + isoform_name
                    exit(1)
                elif (p0 < start_pos[j]):
                    flag = False
                else:
                    j += 1
                    flag = (j < num_exons)

        # Check the junction regions
        for i in range(len(exon_list) - 1):     # the last exon can not have any junction
            p0 = exon_list[i][0]
            p1 = exon_list[i][1] + 1 # the end_pos is not included
            l  = exon_list[i][2]
            region_name_temp = 'P' + str(gene_points_dict[rname][gname][p0]) + ':P' + str(gene_points_dict[rname][gname][p1-1])
            if (region_name_temp not in isoform_exons):   # this exon is not in the isoform
                continue
            # compute start and end point of a region that starts with p1 point
            start = max(p1 - READ_LEN + 1, p0)
            end = p1 - READ_JUNC_MIN_MAP_LEN + min(comp_len_to_junc_gap_forward(i, exon_list, isoform_exons, gene_points_dict[rname][gname]),
                                                   READ_JUNC_MIN_MAP_LEN - 1)
            if (end < start):
                continue    # exon is not long enough to be mapped by a read
            current_point = start
            while (current_point <= end):
                region_name = ""
                if ((current_point == p0) and
                    (p0 != (p1-1))):   # Special case when the exon length is 1 we dont want to repeat the point number
                    region_name += ('P' + str(gene_points_dict[rname][gname][p0]) + '-')
                region_name += ('P' + str(gene_points_dict[rname][gname][p1 - 1]) + '-')
                remain_len = READ_LEN - (p1 - current_point)
                j = i
                while (j < (len(exon_list) -1)):
                    j += 1
                    p0_temp = exon_list[j][0]
                    p1_temp = exon_list[j][1] + 1 # the end_pos is not included
                    l_temp  = exon_list[j][2]
                    region_name_temp = 'P' + str(gene_points_dict[rname][gname][p0_temp]) + ':P' + str(gene_points_dict[rname][gname][p1_temp-1])
                    if (region_name_temp not in isoform_exons):   # this exon is not in the isoform
                        continue
                    if (l_temp >= remain_len):
                        if ((comp_len_to_junc_gap_backward(j, exon_list, isoform_exons, gene_points_dict[rname][gname]) + remain_len ) >= READ_JUNC_MIN_MAP_LEN):
                            region_name += ('P' + str(gene_points_dict[rname][gname][p0_temp]))
                            if ((l_temp == remain_len) and
                                (p0_temp != (p1_temp -1))):   # Special case when the exon length is 1 we dont want to repeat the point number
                                region_name += ('-P' + str(gene_points_dict[rname][gname][p1_temp - 1]))
                            temp_isoform_name = set()
                            temp_isoform_name.add(isoform_name)
                            #temp_isoform_name = {isoform_name}                            
                            if (region_ind.has_key(region_name)):
                                temp_isoform_name = temp_isoform_name.union(region_ind[region_name])
                            else:
                                isoforms_regions_len_dict[rname][gname][region_name] = dict()
                            region_ind[region_name] = temp_isoform_name
                            if (isoforms_regions_len_dict[rname][gname][region_name].has_key(isoform_name)):
                                isoforms_regions_len_dict[rname][gname][region_name][isoform_name] += 1
                            else:
                               isoforms_regions_len_dict[rname][gname][region_name][isoform_name] = 1
                        break   # Found possible region for this current_point
                    else:
                        remain_len -= l_temp
                        region_name += ('P' + str(gene_points_dict[rname][gname][p0_temp]) + '-')
                        if (p0_temp != (p1_temp -1)):   # Special case when the exon length is 1 we dont want to repeat the point number
                            region_name += ('P' + str(gene_points_dict[rname][gname][p1_temp - 1]) + '-')
                current_point += 1
        gene_regions_dict[rname][gname] = region_ind

    sanity_check_isoform_regions_length(gene_exons_dict, gene_regions_dict, gene_isoforms_dict, 
                                        isoforms_regions_len_dict, genes_regions_len_dict)
    
    file_read.close()
    return [gene_exons_dict, gene_points_dict, gene_isoforms_dict, genes_regions_len_dict,
            isoforms_regions_len_dict, gene_regions_dict, gene_isoforms_length_dict]
# #########



# ## Generate gene regions output
# ########################
def generate_regions_output_file(gene_exons_dict, gene_points_dict, gene_isoforms_dict, genes_regions_len_dict,
                                 isoforms_regions_len_dict, gene_regions_dict, gene_isoforms_length_dict,
                                 regions_file_str):
    
    file_write = open(regions_file_str, 'w')
            
    keys1 = sorted(gene_exons_dict.keys())
    for i in keys1:
        keys2 = sorted(gene_exons_dict[i].keys())
        for j in keys2:
            keys3 = sorted(gene_points_dict[i][j].keys())
            keys5 = sorted(gene_isoforms_dict[i][j])

            file_write.write(j.ljust(max(compute_min_ljust(j), 20)))
            file_write.write(str(len(gene_isoforms_dict[i][j])).ljust(10) + i.ljust(10) + '\n')
            for m in keys5:
                file_write.write(m.ljust(max(compute_min_ljust(m), 20)))
            file_write.write('\n')
            for m in keys5:
                file_write.write(str(gene_isoforms_length_dict[i][j + '_' + m]).ljust(max(compute_min_ljust(m), 20)))
            file_write.write('\n')
            for k in keys3:
                file_write.write(('P' + str(gene_points_dict[i][j][k])).ljust(20))
            file_write.write('\n')
            for k in keys3:
                file_write.write(str(k).ljust(20))
            file_write.write('\n')
            keys4 = sorted(gene_regions_dict[i][j])
            for l in keys4:
                file_write.write(l.ljust(max(compute_min_ljust(l), 20)))
            file_write.write('\n')
            for m in keys5:
                for l in keys4:
                    if m in gene_regions_dict[i][j][l]:
                        if (isoforms_regions_len_dict[i][j][l][m] >  0):
                            file_write.write(str(1).ljust(max(compute_min_ljust(l), 20)))
                        else:   # This condition should never be satisfied
                            file_write.write(str(0).ljust(max(compute_min_ljust(l), 20)))
                    else:
                        file_write.write(str(0).ljust(max(compute_min_ljust(l), 20)))
                file_write.write('\n')
            for l in keys4:
                file_write.write(str(genes_regions_len_dict[i][j][l]).ljust(max(compute_min_ljust(l), 20)))
            file_write.write('\n')

    file_write.close()
    
# #########



# ## Main
# #########
def main():
    # Read input parameters
    ref_file_str = sys.argv[1]
    READ_LEN = int(sys.argv[2])
    READ_JUNC_MIN_MAP_LEN = int(sys.argv[3])
    
    # define exon filename
    fields = ref_file_str.split('.')
    exon_file_str = fields[0] + '_exons'
    for i in range(len(fields) - 1):
        exon_file_str += '.' + fields[i + 1]
    
    print 'Generating exon file.' 
    generate_gene_exons(ref_file_str, exon_file_str)
    
    # define isoform filename
    fields = ref_file_str.split('.')
    regions_file_str = fields[0] + '_regions'
    for i in range(len(fields) - 1):
        regions_file_str += '.' + fields[i + 1]
    
    print 'Generating regions file.'
    start_time = datetime.now()
    print 'Generating regions file started at ' + str(start_time)
    
    [gene_exons_dict, gene_points_dict, gene_isoforms_dict, genes_regions_len_dict,
        isoforms_regions_len_dict, gene_regions_dict, gene_isoforms_length_dict] = generate_gene_regions(ref_file_str, exon_file_str, READ_LEN, READ_JUNC_MIN_MAP_LEN)
    generate_regions_output_file(gene_exons_dict, gene_points_dict, gene_isoforms_dict, genes_regions_len_dict,
                                 isoforms_regions_len_dict, gene_regions_dict, gene_isoforms_length_dict, regions_file_str)
    
    end_time = datetime.now()
    print 'Generating regions file ended at   ' + str(end_time)


if __name__ == '__main__':
    main()



