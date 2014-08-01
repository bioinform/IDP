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
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    num_threads = int(sys.argv[3])
    python_path = sys.argv[4]
    penalty_filename = ''
    if (len(sys.argv) > 5):
        penalty_filename = sys.argv[5]
    
    input_file = open(input_filename, 'r' )
    header = input_file.readline()
    input_files = []
    output_filenames = []
    for thread_idx in range(num_threads):
        input_files.append(open(input_filename + '.' + str(thread_idx), 'w'))
        input_files[-1].write(header)

        
    thread_idx = 0
    while True:
        line = input_file.readline()
        if line == "": 
            break
        num_isoforms = int(line.split()[1])
        input_files[thread_idx].write(line)
        
        for i in range(6):
            input_files[thread_idx].write(input_file.readline())
        for i in range(num_isoforms):
            input_files[thread_idx].write(input_file.readline())
        for i in range(2):
            input_files[thread_idx].write(input_file.readline())
    
        thread_idx = (thread_idx + 1) % num_threads
        
    for thread_idx in range(num_threads):
        input_files[thread_idx].close()
    input_file.close()
    
    ##############################
    threads_list = []
    for thread_idx in range(num_threads):
        cmd = python_path + " " + bin_path + 'MLE_regions.py ' + input_filename + '.' + str(thread_idx) + ' ' + output_filename + '.' + str(thread_idx)
        if (penalty_filename != ''):
            cmd += ' ' + penalty_filename
        print cmd
        threads_list.append( threading.Thread(target=os.system, args=(cmd,)) )
        threads_list[thread_idx].start()

    for thread in threads_list:
        thread.join()
        
    cat_cmnd = 'cat '
    rm_cmnd = "rm "
    for thread_idx in range(num_threads):
        cat_cmnd += output_filename + '.' + str(thread_idx) + " "
        rm_cmnd += output_filename + '.' + str(thread_idx)  + " " + input_filename + '.' + str(thread_idx) + " "
        
    cat_cmnd += ' > ' + output_filename
    os.system(cat_cmnd)
    os.system(rm_cmnd)
    
        

if __name__ == '__main__':
    main()
