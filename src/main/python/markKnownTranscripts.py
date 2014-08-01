#!/usr/bin/python
import sys

### 
##########
def compute_min_ljust(str_value):
    return max(20, (len(str_value)/10 + 1) * 10)


###
###########
def parse_known_transcripts(transcript_filename):
    
    transc_dict = {}
        
    transcript_file = open(transcript_filename, 'r')
    for line in transcript_file:
        fields = line.split()
        if not transc_dict.has_key(fields[0]):
            transc_dict[fields[0]] = set()
        transc_dict[fields[0]].add(fields[1])
    
    transcript_file.close()
    return transc_dict

###
#########
def mark_known_transcripts(input_filename, output_filename,
                           transc_dict):
    
    input_file = open(input_filename, 'r')
    output_file = open(output_filename, 'w')
    
    output_file.write(input_file.readline())
    line_index = 0;
    for line in input_file:
        fields = line.split()
        output_file.write(line)
        if (line_index == 0):
            gname = fields[0]
            num_isoforms = int(fields[1])
            rname = fields[2]
        if (line_index == 1):
            for i in range(num_isoforms):
                if (transc_dict.has_key(rname) and (fields[i] in transc_dict[rname])):
                    output_file.write(str(1).ljust(compute_min_ljust(fields[i])))
                else:
                    output_file.write(str(0).ljust(compute_min_ljust(fields[i])))
            output_file.write('\n')
        line_index += 1
        if (line_index == (8 + num_isoforms)):
            line_index = 0
            
    input_file.close()
    output_file.close()
    
    
### Main
##########
def main():
    # Read input parameters

    input_filename = sys.argv[1]
    transcript_filename = sys.argv[2]
    output_filename = sys.argv[3]
    
    transc_dict = parse_known_transcripts(transcript_filename)
    mark_known_transcripts(input_filename, output_filename, transc_dict)

if __name__ == '__main__':
    main()




