#!/usr/bin/python

import sys
import os
import threading

def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    return path, filename

#Main**************************************************************************#
def main():
    # Read input parameters
    bin_path,command = GetPathAndName(sys.argv[0])
    regions_filename = sys.argv[1]
    reads_filename = sys.argv[2]
    num_threads = int(sys.argv[3])
    python_path = sys.argv[4]
    read_len = sys.argv[5]
    min_junction_overlap_len = sys.argv[6]
    output_filename = 'refSeq_MLE_input.txt'
    
    reads_file = open(reads_filename, 'r' )
    reads_files = []
    for thread_idx in range(num_threads):
        reads_files.append(open(reads_filename + '.' + str(thread_idx), 'w'))

        
    thread_idx = 0
    for line in reads_file:
        reads_files[thread_idx].write(line)
        thread_idx = (thread_idx + 1) % num_threads
        
    for thread_idx in range(num_threads):
        reads_files[thread_idx].close()
    reads_file.close()
    
    ##############################

    threads_list = []
    for thread_idx in range(num_threads):
        cmd = (python_path + " " + bin_path + 'parseSAM.py ' + regions_filename + ' ' + reads_filename + '.' + str(thread_idx) + 
               ' ' + output_filename + '.' + str(thread_idx) + ' ' + read_len + ' ' + min_junction_overlap_len)
        print cmd
        threads_list.append( threading.Thread(target=os.system, args=(cmd,)) )
        threads_list[thread_idx].start()

    for thread in threads_list:
        thread.join()
    
    output_file = open(output_filename, 'w')
    output_files = []
    
    header = 0
    for thread_idx in range(num_threads):
        output_files.append(open(output_filename + '.' + str(thread_idx), 'r'))
        header += int(output_files[thread_idx].readline())
    output_file.write(str(header) + '\n')
    
    genes_str_map = {}
    genes_reads_count_map = {}
    
    for thread_idx in range(num_threads):
        while True:
            line = output_files[thread_idx].readline()
            if (line == ''):
                break

            if not genes_str_map.has_key(line):
                lines = line 
                isoforms_line = output_files[thread_idx].readline()
                lines += isoforms_line 
                for i in range(4):
                    lines += output_files[thread_idx].readline()
                for i in range(len(isoforms_line.split())):
                    lines += output_files[thread_idx].readline()
                lines += output_files[thread_idx].readline()
                genes_reads_count_map[line] = [int(i) for i in output_files[thread_idx].readline().split()]

                genes_str_map[line] = lines
            else:
                isoforms_line = output_files[thread_idx].readline()
                for i in range(4):
                    output_files[thread_idx].readline()
                for i in range(len(isoforms_line.split())):
                    output_files[thread_idx].readline()
                output_files[thread_idx].readline()
                reads_count_line = [int(i) for i in output_files[thread_idx].readline().split()]
                for i in range(len(reads_count_line)):
                    genes_reads_count_map[line][i] += reads_count_line[i]

                
    genes = sorted(genes_str_map.keys())
    for gene in genes:
        output_file.write(genes_str_map[gene])
        for i in range(len(genes_reads_count_map[gene])):
            output_file.write(str(genes_reads_count_map[gene][i]).ljust(20))
        output_file.write('\n')
        
    
    
    for thread_idx in range(num_threads):
        output_files[thread_idx].close()
        rm_cmnd = "rm " + output_filename + '.' + str(thread_idx) + ' ' + reads_filename + '.' + str(thread_idx)
        os.system(rm_cmnd)
    output_file.close()
    
    
if __name__ == '__main__':
    main()
