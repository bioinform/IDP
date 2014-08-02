#!/usr/bin/python
import sys
import struct
import os
from numpy import *
import datetime
from re import *
from copy import *
import threading
import string
from binaidp import log_command
import inspect # JW for debugging

################################################################################


######### GetPathAndName ##########
# Break a path+filename into its components
#
# Pre:  Takes a filename
# Post: Return the path and filename strings
# Modifies:  None
###################################
def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    return path, filename


######### Readcfgfile ##########
# Read/Parse in the configuration file
#
# Pre:  Takes a filename
# Post: Return a dictionary of arguments and values
# Modifies:  Standard Out
###################################
def Readcfgfile(cfg_filename):
    results = {}
    cfg = open(cfg_filename,'r')
    for line in cfg:
        line = line.strip()
        if line=='':
            continue
        if not line[0]=='#':
            ls = line.split('=')
            if len(ls)>2:
                print 'warning: too many = in cfg file'
            results[ls[0].strip()] = ls[1].strip()
            print [ls[0].strip(), ls[1].strip()]

    cfg.close()
    return results

######### print_run ##########
# Print and run a command both
#
# Pre:  Takes a command as a string
# Post: 
# Modifies:  Standard Out, Calls os.system()
###################################
def print_run(cmd, ignorefail=False):
    print cmd
    print ""
    log_command(cmd, ignorefail)

######### folder_absolute_right_slash ##########
# Make sure the path for a folder has right slash and is an absolute path
#
# Pre:  Takes a path as a string, Uses os.path.abspath
# Post: Outputs a path as a string
# Modifies:  None
###################################
def folder_absolute_right_slash(path):
  path = os.path.abspath(path)
  if path[-1]!='/':
    path=path+'/'
  return path

######### debug_stop ##########
# Exit and print a message stating the location of the stop
#
# Pre:  String, Requires the inspect module
# Post: 
# Modifies:  Standard Out, calls sys.exit()
###################################
def debug_stop(var):
    print "Debug exit at line " + str(inspect.currentframe().f_back.f_lineno) + "."
    if var:
      print "Debug message: " + var
    sys.exit()

######### debug_report ##########
# Print a message stating the execution of a line
#
# Pre:  String, Requires the inspect module
# Post: 
# Modifies:  Standard Out
###################################
def debug_report(var):
    print "Debug passing line " + str(inspect.currentframe().f_back.f_lineno) + "."
    if var:
      print "Debug message: " + var



################################################################################

# Check for required user inputs
if len(sys.argv) >= 3:
    run_pathfilename =  sys.argv[0]
    cfg_filename =  sys.argv[1]
    Istep = int(sys.argv[2])
else:
    print "usage: python runIDP.py run.cfg 3"
    print "or ./runIDP.py run.cfg 3"
    sys.exit(1)

################################################################################
# Initialize variables from configuration
python_path = "/usr/bin/python"
LR_gpd_pathfilename = ""
LR_psl_pathfilename = ""
SR_jun_pathfilename = ""
SR_sam_pathfilename = ""
CAGE_data_filename = ""
ref_annotation_pathfilename = ""
allref_annotation_pathfilename = ""
Njun_limit = "20"
L_exon_limit = "4000"
Ijunction_cover = "1"
temp_foldername = ""
output_foldername = ""
detected_exp_len = ""
Bfile_Npt = "500"
Nbin = "5"
I_refjun_isoformconstruction =  "1"
I_ref5end_isoformconstruction = "1"
I_ref3end_isoformconstruction = "1"

three_primer = ""
five_primer = ""
exon_construction_junction_span = "1"
Niso = "50"
Npt = "500"
Nbin = "5"
FPR = "0.05"

read_length = 100
L_min_intron = 68
min_junction_overlap_len = 10
psl_type = "0"
################################################################################
cfg_dt = Readcfgfile(cfg_filename)
for key in cfg_dt:   # Assign variables from configuration file
    if key == "Nthread":
        Nthread = int(cfg_dt[key])
    elif key == "python_path":
        python_path = cfg_dt[key]
    elif key == "LR_gpd_pathfilename":
        LR_gpd_pathfilename = cfg_dt[key]
    elif key == "SR_jun_pathfilename":
        SR_jun_pathfilename = os.path.abspath(cfg_dt[key])
    elif key == "SR_sam_pathfilename":
        SR_sam_pathfilename = os.path.abspath(cfg_dt[key])
    elif key == "ref_annotation_pathfilename":
        ref_gpd_pathfilename = os.path.abspath(cfg_dt[key])
    elif key == "allref_annotation_pathfilename":
        allref_gpd_pathfilename = os.path.abspath(cfg_dt[key])
    elif key == "min_junction_overlap_len":
        min_junction_overlap_len = int(cfg_dt[key])
    elif key == "read_length":
        read_length = int(cfg_dt[key])
        
    elif key == "Njun_limit":
        Njun_limit = cfg_dt[key]
    elif key == "L_exon_limit":
        L_exon_limit = cfg_dt[key]
    elif key == "Niso_limit":
        Niso = cfg_dt[key]
    elif key == "Ijunction_cover":
        Ijunction_cover = cfg_dt[key]
    elif key == "L_min_intron":
        L_min_intron = int(cfg_dt[key])

    elif key == "detected_exp_len":
        detected_exp_len_pathfilename = os.path.abspath(cfg_dt[key])
    elif key == "Bfile_Npt":
        Npt = cfg_dt[key]
    elif key == "Bfile_Nbin":
        Nbin = cfg_dt[key]
    elif key == "exon_construction_junction_span":
        exon_construction_junction_span =cfg_dt[key]

    elif key == "temp_foldername":
        temp_foldername = cfg_dt[key]
    elif key == "output_foldername":
        output_foldername = cfg_dt[key]
    elif key == "I_refjun_isoformconstruction":
        I_refjun_isoformconstruction = cfg_dt[key]
    elif key == "I_ref5end_isoformconstruction":
        I_ref5end_isoformconstruction = cfg_dt[key]
    elif key == "I_ref3end_isoformconstruction":
        I_ref3end_isoformconstruction = cfg_dt[key]

    elif key == "CAGE_data_filename":
        CAGE_data_filename = cfg_dt[key]

    elif key == "LR_pathfilename":
        LR_pathfilename = cfg_dt[key]
    elif key == "genome_pathfilename":
        genome_pathfilename = cfg_dt[key]
    elif key == "LR_psl_pathfilename":
        LR_psl_pathfilename = cfg_dt[key]
    elif key == "psl_type":
        psl_type = cfg_dt[key]

    elif key == "three_primer":
        three_primer = cfg_dt[key]
    elif key == "five_primer":
        five_primer = cfg_dt[key]

    elif key == "FPR":
        FPR = cfg_dt[key]

################################################################################
# Create folders in the paths supplied in the configuration file
print_run('mkdir -pv ' + temp_foldername)
print_run('mkdir -pv ' + output_foldername)

# Identify the folder containing the binaries
bin_foldername, run_filename = GetPathAndName(run_pathfilename)

temp_foldername = folder_absolute_right_slash(temp_foldername)
output_foldername = folder_absolute_right_slash(output_foldername)
bin_foldername = folder_absolute_right_slash(bin_foldername)
python_bin_foldername = python_path + " " + bin_foldername
begin_dir = os.getcwd()

##check input from cfg file######################################################

## CAGE data ##
if I_ref5end_isoformconstruction != "1" and CAGE_data_filename == "":
    print "Error: I_ref5end_isoformconstruction != 1","no CAGE data"
    sys.exit(1)

## long reads ##
# Priority gpd > psl > fa/blat
I_LR_step = 0
if LR_gpd_pathfilename != "":
    print "use the BLAT alignment of long reads (gpd format), " + LR_gpd_pathfilename + " as long read input"
    I_LR_step = 0
    print_run("cp " + LR_gpd_pathfilename +  " " + temp_foldername + "LR.gpd")

elif LR_psl_pathfilename != "":
    print "use the BLAT alignment of long reads (psl format), " + LR_psl_pathfilename + " as long read input"
    I_LR_step = 1

    if psl_type == "0":
        blat_best_cmd = python_bin_foldername + "blat_best.py " + LR_psl_pathfilename + " 5 > " + temp_foldername + "LR.bestpsl"
        print_run(blat_best_cmd)
    else:
        blat_best_cmd = "cp " + LR_psl_pathfilename + " " + temp_foldername + "LR.bestpsl"
        print_run(blat_best_cmd)

    change_psl_cmd = python_bin_foldername + "change_psl_4digit.py " + temp_foldername + "LR.bestpsl > " + temp_foldername + "newname4_LR.bestpsl"
    print_run(change_psl_cmd)

    psl2genephed_cmd = python_bin_foldername + "psl2genephed.py " + temp_foldername + "newname4_LR.bestpsl 0 " + str(L_min_intron) + " " + temp_foldername + "LR.gpd"
    print_run(psl2genephed_cmd)
	
elif LR_pathfilename != "" and genome_pathfilename != "":
    print "use raw sequence of the long reads (FASTA format), " + LR_pathfilename + " as long read input; and the reference genome is " + genome_pathfilename
    #This sectino will make a large number of files in the temp folder (~156)
    I_LR_step = 2
    if three_primer != "" and five_primer != "":
        compress_cmd = python_bin_foldername + "compressFASTA.py " + LR_pathfilename +  " " + temp_foldername + "LR."
        print_run(compress_cmd)

        print_run( "echo \">three_primer\" > " + temp_foldername + "three.fa" )
        print_run( "echo \"" + three_primer + "\" >> " + temp_foldername + "three.fa" )

        seqmap_three_cmd = bin_foldername + "seqmap 2 " + temp_foldername + "three.fa "  + temp_foldername + "LR.cps " + temp_foldername + "LR.cps.3.out /allow_insdel:1  /output_alignment /output_all_matches"
        print_run(seqmap_three_cmd)

        print_run( "echo \">five_primer\" > " + temp_foldername + "five.fa" )
        print_run( "echo \"" + five_primer + "\" >> " + temp_foldername + "five.fa" )
	
        seqmap_five_cmd = bin_foldername + "seqmap 2 " + temp_foldername + "five.fa "  + temp_foldername + "LR.cps " + temp_foldername + "LR.cps.5.out /allow_insdel:1  /output_alignment /output_all_matches"
        print_run(seqmap_five_cmd)

        removeAdapterPolyA_cmd = python_bin_foldername + "removeAdapterPolyA.py " + temp_foldername + "LR.cps " + temp_foldername + "LR.idx " + temp_foldername + "LR.cps.3.out " + temp_foldername + "LR.cps.5.out " + temp_foldername + "LR_notailspolyA.fa"
        print_run(removeAdapterPolyA_cmd)

        blat_cmd = python_bin_foldername + "blat_threading.py " + python_path + " " + str(Nthread) + " -t=DNA -q=DNA " + genome_pathfilename + " " + temp_foldername + "LR_notailspolyA.fa " + temp_foldername + "LR_notailspolyA.fa.bestpsl" #JWDEBUG added a python path argument to fix blat_threading.py's call to best_blat.py
        print_run(blat_cmd)
		
        change_psl_cmd = python_bin_foldername + "change_psl_polyA3end_4digit.py " + temp_foldername + "LR_notailspolyA.fa.bestpsl " + temp_foldername + "LR_notailspolyA.fa.3 > " + temp_foldername + "newname4_LR.bestpsl"
        print_run(change_psl_cmd)

        psl2genephed_cmd = python_bin_foldername + "psl2genephed.py " + temp_foldername + "newname4_LR.bestpsl 0 " + str(L_min_intron) + " " + temp_foldername + "LR.gpd"
        print_run(psl2genephed_cmd)
	
    else:
        blat_cmd = python_bin_foldername + "blat_threading.py " + python_path + " " + str(Nthread) + " -t=DNA -q=DNA " + genome_pathfilename + " " + LR_pathfilename + " " + temp_foldername + "LR.fa.psl" #JWDEBUG also did it to this line
        print_run(blat_cmd)
		
        change_psl_cmd = python_bin_foldername + "change_psl_4digit.py " + temp_foldername + "LR.fa.psl > " + temp_foldername + "newname4_LR.bestpsl"
        print_run(change_psl_cmd)

        psl2genephed_cmd = python_bin_foldername + "psl2genephed.py " + temp_foldername + "newname4_LR.bestpsl 0 " + str(L_min_intron) + " " + temp_foldername + "LR.gpd"
        print_run(psl2genephed_cmd)

else:
    print "no long read data (alignment or raw sequences)"
    sys.exit(1)

## short reads ##
if SR_sam_pathfilename != "" and SR_jun_pathfilename != "":
    print "use the short read alignment (SAM format) " + SR_sam_pathfilename + " and junctions (bed format) " + SR_jun_pathfilename + " as short read input"
else:
    print "Error: no short read alignment data (SR_sam_pathfilename) or junction detection (SR_jun_pathfilename)"
    print "Run SpliceMap (http://www.stanford.edu/group/wonglab/SpliceMap/) or Tophat to generate the read alignment data and junction detection"
    sys.exit(1)

################################################################################

print_run("cp " + SR_jun_pathfilename + " " + temp_foldername + "SR.bed")
print_run("cp " + ref_gpd_pathfilename + " " + temp_foldername + "ref.gpd")

#################################################################################

## detected_exp_len ##
I_sam_exist = 0
try:
    detected_exp_len_pathfilename
except NameError:
    print "Warning: There is no " + detected_exp_len + "data." 
    print "Here, we calculate detection rate from long reads data and short read alignment" + SR_sam_pathfilename

    os.chdir(temp_foldername)
    
    ## abundance estimation of annotated transcripts ##

    parseRef_cmd = python_bin_foldername + "parseRef.py " + "ref.gpd " + str(read_length) + " " + str(min_junction_overlap_len)
    print_run(parseRef_cmd)

    if SR_sam_pathfilename.endswith('.bam'):
        print_run("samtools view -h -o SR.sam " + SR_sam_pathfilename)
    else:
        print_run("cp " + SR_sam_pathfilename + " SR.sam")
    I_sam_exist = 1

    parseSAM_cmd = python_bin_foldername + "parseSAM_MT.py " + "ref_regions.gpd " + "SR.sam " + str(Nthread) + " " + python_path + " " +  + str(read_length) + " " + str(min_junction_overlap_len) + " > parseSAM_MT0.log" 
    print_run(parseSAM_cmd)
    print_run("awk \'{print $3\"\\t\"$2}\' " + "ref.gpd > " + "positive_candidate_list0")

    print_run("mv refSeq_MLE_input.txt refSeq_MLE_input0.txt")

    markknownTranscripts_cmd = python_bin_foldername + "markKnownTranscripts.py " + "refSeq_MLE_input0.txt " + "positive_candidate_list0 " + "refSeq_MLE_input_marked0.txt" 
    print_run(markknownTranscripts_cmd)

    MLE_cmd = python_bin_foldername + "MLE_MT.py " + "refSeq_MLE_input_marked0.txt " + "refSeq_MLE_output0.txt " + str(Nthread) + " " + python_path
    print_run(MLE_cmd)

    os.chdir(begin_dir)

    ####################################################


    novel_genephed_cmd = python_bin_foldername + "novel_genephed.py " + temp_foldername + "ref.gpd" + " " + temp_foldername + "LR.gpd " + temp_foldername + "novel_LR.gpd > " + temp_foldername + "known_LR.gpd"
    print_run(novel_genephed_cmd)
   
    ## abundance estimation of annotated transcripts ##

    reformat_cmd = python_bin_foldername + "reformat.py " + temp_foldername + "refSeq_MLE_output0.txt > " + temp_foldername + "refSeq_MLE_output0.txt_"
    print_run(reformat_cmd)

    maketab_cmd = python_bin_foldername + "MLEout2tab.py " + temp_foldername + "refSeq_MLE_output0.txt_ > " + temp_foldername + "refSeq_MLE_output0.tab"
    print_run(maketab_cmd)

    exp_len_I_cmd = python_bin_foldername + "exp_len.py " + temp_foldername + "refSeq_MLE_output0.tab " + temp_foldername + "known_LR.gpd_ref.gpd > " + temp_foldername + "known_LR.gpd_ref.gpd_exp_len"
    print_run(exp_len_I_cmd)

else:
    print_run("cp " + detected_exp_len_pathfilename + " " + temp_foldername + "known_LR.gpd_ref.gpd_exp_len")

#############################################################################################################################################################

if Istep == 1 or Istep == 0:
    combine_jun_cmd = python_bin_foldername + "gpd2jun.py " + temp_foldername + "ref.gpd" + " > " + temp_foldername + "ref.gpd.bed"
    print_run(combine_jun_cmd)

    print_run("cat " + temp_foldername + "ref.gpd.bed " + temp_foldername + "SR.bed" + " > " + temp_foldername + "ref_experiment_jun.bed ")

    juncover_cmd = python_bin_foldername + "juncover_gpd.py " + temp_foldername + "ref_experiment_jun.bed " + temp_foldername + "LR.gpd " + temp_foldername + "filterout_LR_jun.bed > " + temp_foldername + "junfil_LR.gpd"
    print_run(juncover_cmd)

################################################################################

    exon_construction_cmd = python_bin_foldername + "exon_construction.py " + temp_foldername + "ref.gpd" + " " + temp_foldername + "SR.bed " + exon_construction_junction_span + " > " + temp_foldername + "novel_long_junctions.bed"
    print_run(exon_construction_cmd)

    compatible_cmd = python_bin_foldername + "compatible_LR_2_polyA3end.py " + temp_foldername + "junfil_LR.gpd " + temp_foldername + "SR.bed_ref.gpd.exon > " + temp_foldername + "junfil_compatible_LR_polyA3end.log"
    print_run(compatible_cmd)

    print_run( "grep -v end " + temp_foldername + "junfil_compatible_LR_polyA3end.log | grep -v termin > " + temp_foldername + "junfil_compatible_LR_polyA3end.gpd" )

    jun_isoformconstruction = temp_foldername + "SR.bed"
#    if I_refjun_isoformconstruction == "1":
#         jun_isoformconstruction  = temp_foldername + "ref_experiment_jun.bed"

#    isoform_construction_cmd = bin_foldername + "isoform_construction_polyA3end_5cap.py " + temp_foldername + "ref.gpd" + " " + jun_isoformconstruction + " " + temp_foldername + "SR.bed" + "_ref.gpd.exon " + CAGE_data_filename + " " + temp_foldername + "junfil_compatible_LR_polyA3end.gpd" + " " + I_ref5end_isoformconstruction + " " + I_ref3end_isoformconstruction + " " + Njun_limit + " " + L_exon_limit + " > " + temp_foldername + "isoform_construction.log"
#    print_run(isoform_construction_cmd)

    ##############################CAGE

#    python2.6 mapEncodeTSStoRegions.py encodeTssHmm.bedRnaElements SR.bed_ref.gpd.exon
    processed_CAGE_filename = "kinfai"
    if CAGE_data_filename != "":
        curr_dir = os.getcwd()
        os.chdir(temp_foldername)
        experimental_5end_cmd = python_bin_foldername + "mapEncodeTSStoRegions.py " + CAGE_data_filename + " " + "SR.bed_ref.gpd.exon"
        print_run(experimental_5end_cmd)
        os.chdir(curr_dir)
        print_run(python_bin_foldername + "reformat.py " + temp_foldername + "encodeTSS_mapped_regions.txt > " + temp_foldername + "encodeTSS_mapped_regions.txt_")
        processed_CAGE_filename = temp_foldername + "encodeTSS_mapped_regions.txt_"

    ##############################

    temp_LR_gpd = open(temp_foldername + "junfil_compatible_LR_polyA3end.gpd",'r')
    temp_LR_gpd_NR = 0
    for line in temp_LR_gpd:
        temp_LR_gpd_NR+=1
    temp_LR_gpd.close()
    
    if (temp_LR_gpd_NR == 0):
        print "Nothing found in the junfil_compatible_LR_polyA3end.gpd file!";
        sys.exit(1);

    Nsplitline = temp_LR_gpd_NR/Nthread
    print Nsplitline

    ext_ls=[]
    for i in range(Nthread):
        ext_ls.append( '.' + string.lowercase[i/26] + string.lowercase[i%26] )


    split_LR_gpd_cmd = "split -l " + str(Nsplitline) + " " + temp_foldername + "junfil_compatible_LR_polyA3end.gpd" + " " + temp_foldername + "junfil_compatible_LR_polyA3end.gpd" +"."
    print_run(split_LR_gpd_cmd)

    ##############################

    i = 0
    T_isoform_construction_ls = []
    for ext in ext_ls:
        isoform_construction_cmd = python_bin_foldername + "isoform_construction_polyA3end_5cap.py " + temp_foldername + "ref.gpd" + " " + jun_isoformconstruction + " " + temp_foldername + "SR.bed" + "_ref.gpd.exon " + processed_CAGE_filename + " " + temp_foldername + "junfil_compatible_LR_polyA3end.gpd" + ext + " " + I_ref5end_isoformconstruction + " " + I_ref3end_isoformconstruction + " " + Njun_limit + " " + L_exon_limit + " " + temp_foldername + "junfil_compatible_LR_polyA3end.gpd.out" + ext + " > " + temp_foldername + "isoform_construction" + ext + ".log"
        print isoform_construction_cmd
        T_isoform_construction_ls.append( threading.Thread(target=os.system, args=(isoform_construction_cmd,)) )
        T_isoform_construction_ls[i].start()
        i+=1
    for T in T_isoform_construction_ls:
        T.join()

    merge_isoform_construction_cmd = python_bin_foldername + "merge_isoform_construction_polyA3end_5cap.py " + temp_foldername + "isoform_construction.gpd "
    for ext in ext_ls:
        merge_isoform_construction_cmd = merge_isoform_construction_cmd + " " + temp_foldername + "junfil_compatible_LR_polyA3end.gpd.out" + ext
    print_run(merge_isoform_construction_cmd)

#    print_run( "grep -v warning " + temp_foldername + "isoform_construction.log | grep -v ending > " + temp_foldername + "isoform_construction.log2")

################################################################################
    print_run("cut -f 1 " + temp_foldername + "isoform_construction.gpd > " + temp_foldername + "isoform_construction.gpd.gene")
    print_run("sort " + temp_foldername + "isoform_construction.gpd.gene > " + temp_foldername + "isoform_construction.gpd.gene_sorted" )
    print_run("uniq -c " + temp_foldername + "isoform_construction.gpd.gene_sorted > " + temp_foldername + "isoform_construction.gpd.gene_sorted.count")
    print_run(python_bin_foldername + "reformat.py " + temp_foldername + "isoform_construction.gpd.gene_sorted.count > " + temp_foldername + "isoform_construction.gpd.gene_sorted.count_")
    print_run( "awk \'{ if($1<=" + Niso + ") print $2}\' " + temp_foldername + "isoform_construction.gpd.gene_sorted.count_ > " + temp_foldername + "isoform_construction.gpd." + "Niso" + Niso  + ".gene")
    print_run( python_bin_foldername + "selectrow.py " + temp_foldername + "isoform_construction.gpd" + " 1 " + temp_foldername + "isoform_construction.gpd." + "Niso" + Niso  + ".gene" + " " + temp_foldername + "isoform_construction." + "Niso" + Niso  + ".gpd > " + temp_foldername + "isoform_construction.gpd." + "Niso" + Niso  + ".gpd_selectcount")

################################################################################

#############################################################################################################################################################
if Istep == 2 or Istep == 0:
    os.chdir(temp_foldername)

########################
#    print_run( "grep chr19 " + "isoform_construction." + "Niso" + Niso + ".gpd" + " > " + "chr19_isoform_construction." + "Niso" + Niso + ".gpd" )
#    print_run( "mv isoform_construction." + "Niso" + Niso + ".gpd" + " " + "all_isoform_construction." + "Niso" + Niso + ".gpd")    
#    print_run( "mv chr19_isoform_construction." + "Niso" + Niso + ".gpd" + " " + "isoform_construction." + "Niso" + Niso + ".gpd")    
########################

    parseRef_cmd = python_bin_foldername + "parseRef.py " + "isoform_construction." + "Niso" + Niso + ".gpd " + str(read_length) + " " + str(min_junction_overlap_len)
    print_run(parseRef_cmd)

    if I_sam_exist == 0:
        if SR_sam_pathfilename.endswith('.bam'):
            print_run("samtools view -h -o SR.sam " + SR_sam_pathfilename)
        else:
            print_run("cp " + SR_sam_pathfilename + " SR.sam")

    parseSAM_cmd = python_bin_foldername + "parseSAM_MT.py " + "isoform_construction_regions." + "Niso" + Niso  + ".gpd " + "SR.sam " + str(Nthread) + " " + python_path + " " + str(read_length) + " " + str(min_junction_overlap_len) + " > parseSAM_MT.log"
    print_run(parseSAM_cmd)
    novel_genephed_cmd = python_bin_foldername + "novel_genephed.py " + "ref.gpd" + " " + "LR.gpd " + "novel_LR.gpd > " + "known_LR.gpd"
    print_run(novel_genephed_cmd)

    positive_candidate_cmd = python_bin_foldername + "novel_genephed.py " + "isoform_construction.gpd " + "known_LR.gpd_ref.gpd " + "negative_known_LR.gpd_ref.gpd_isoform_construction.gpd > " + "nonnegative_known_LR.gpd_ref.gpd_isoform_construction.gpd "
    print_run(positive_candidate_cmd)

    print_run("awk '{print $3\"\\t\"$2}' " + "known_"+ "known_LR.gpd_ref.gpd_" + "isoform_construction.gpd > " + "positive_candidate_list")

    markknownTranscripts_cmd = python_bin_foldername + "markKnownTranscripts.py " + "refSeq_MLE_input.txt " + "positive_candidate_list " + "refSeq_MLE_input_marked.txt" 
    print_run(markknownTranscripts_cmd)

    exp_len_I_cmd = python_bin_foldername + "exp_len_on_refbigtable.py " + "known_LR.gpd_ref.gpd_exp_len" + " " + "known_LR.gpd_ref.gpd  > " + "known_LR.gpd_ref.gpd_exp_len_I"
    print_run(exp_len_I_cmd)

    Bfile_cmd = python_bin_foldername + "Bfile.py " + "known_LR.gpd_ref.gpd_exp_len_I " + "isoform_construction." + "Niso" + Niso + ".gpd " + Npt + " " + Nbin + " > " + "Bfile"
    print_run(Bfile_cmd)

    MLE_cmd = python_bin_foldername + "MLE_MT.py " + "refSeq_MLE_input_marked.txt " + "refSeq_MLE_output.txt " + str(Nthread) + " " + python_path + " Bfile "

    print_run(MLE_cmd)
    reformat_cmd = python_bin_foldername + "reformat.py " + "refSeq_MLE_output.txt > " + "refSeq_MLE_output.txt_"
    print_run(reformat_cmd)

    maketab_cmd = python_bin_foldername + "MLEout2tab.py " + "refSeq_MLE_output.txt_ > " + "refSeq_MLE_output.tab"
    print_run(maketab_cmd)

    print_run("mv refSeq_MLE_input.txt refSeq_MLE_input.txt_")

#############################################################################################################################################################

    print_run("cut -f 2 " + "positive_candidate_list > " + "positive_candidate_ID_list")

    selectrow_cmd = python_bin_foldername + "selectrow.py " + "refSeq_MLE_output.tab 1 " + "positive_candidate_ID_list" + " " + "positive_refSeq_MLE_output.tab > positive_candidate_ID_list_selectcount" #JWDEBUG modified this to use the user defined python path
    print_run(selectrow_cmd)

    consective_junlk_filter_cmd = python_bin_foldername + "consective_junlk_filter.py " + allref_gpd_pathfilename + " " + "isoform_construction.gpd " + "isoform_construction.gpd" + " > single-junction_transcript"
    print_run(consective_junlk_filter_cmd)

    print_run( "cut -f 2 " + "negative_isoform_construction.gpd > " + "negative_candidate_ID_list" )

    selectrow_cmd = python_bin_foldername + "selectrow.py " + "refSeq_MLE_output.tab 1 " + "negative_candidate_ID_list" + " " + "negative_refSeq_MLE_output.tab > negative_candidate_ID_list_selectcount"
    print_run(selectrow_cmd)

    addexp2bed_cmd = python_bin_foldername + "addexp2bed.py " + "refSeq_MLE_output.txt_ " + "isoform_construction.gpd isoform_construction_expcol.bed B > noexp_gene_transcript_list"
    print_run(addexp2bed_cmd)

    select_ROC_cmd = python_bin_foldername + "select_FPR.py positive_refSeq_MLE_output.tab negative_refSeq_MLE_output.tab " + str(FPR) + " refSeq_MLE_output.tab refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".tab"
    print_run(select_ROC_cmd)
    print_run("cut -f 1 refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".tab > refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".ID")
    print_run(python_bin_foldername + "selectrow.py " + "isoform_construction." + "Niso" + Niso + ".gpd" +  " 2 " + "refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".ID " + "refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".gpd > " + "isoform_construction." + "Niso" + Niso + ".gpd_selectcount")

    ## Merge Detection and Prediction ##

    assign_genenametoPrediction_junusage_cmd = python_bin_foldername + "assign_genenametoPrediction_junusage.py ref.gpd " + "refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".gpd " + "refname_refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".gpd"
    print_run(assign_genenametoPrediction_junusage_cmd)

    mergeDPgpd_cmd = python_bin_foldername + "mergeDPgpd.py " + "refname_refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".gpd " + "refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".gpd " + "Notref_refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".gpd "
    print_run(mergeDPgpd_cmd)

    print_run("awk \'{if($4!=\"-\")print $4}\' " + "refname_refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".gpd > " + "refname_refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".ID")

    print_run("cut -f 2 known_LR.gpd_ref.gpd > known_LR.gpd_ref.gpd_ID")

    selectrow_cmd = python_bin_foldername + "selectrow.py ref.gpd 2 known_LR.gpd_ref.gpd_ID known_LR.gpd_ref.gpd_ID_ref.gpd > " + "known_LR.gpd_ref.gpd_ID_selectcount" 
    print_run(selectrow_cmd)

    print_run("cut -f 2 known_LR.gpd_ref.gpd_ID_ref.gpd > known_LR.gpd_ref.gpd_ID_ref.gpd_ID")

    union_cmd = python_bin_foldername + "union.py known_LR.gpd_ref.gpd_ID_ref.gpd_ID " + "refname_refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".ID" + " | cut -f 1 > ref_DP_ID"
    print_run(union_cmd)

    selectrow_cmd = python_bin_foldername + "selectrow.py ref.gpd 2 ref_DP_ID ref_DP.gpd > ref_DP_ID_selectcount"
    print_run(selectrow_cmd)

    print_run("cat ref_DP.gpd Notref_refSeq_MLE_output_FPR" +str(FPR).replace(".","") + ".gpd > DP" + str(FPR).replace(".","") + ".gpd")

    ## abundance estimation of ouput transcripts ##

    parseRef_cmd = python_bin_foldername + "parseRef.py " + "DP" +  str(FPR).replace(".","") + ".gpd " + str(read_length) + " " + str(min_junction_overlap_len)
    print_run(parseRef_cmd)

    parseSAM_cmd = python_bin_foldername + "parseSAM_MT.py " + "DP" +  str(FPR).replace(".","") + "_regions.gpd " + "SR.sam " + str(Nthread) + " " + python_path + " "  + str(read_length) + " " + str(min_junction_overlap_len) + " > parseSAM_MT1.log"
    print_run(parseSAM_cmd)

    print_run("awk \'{print $3\"\\t\"$2}\' " + "DP" +  str(FPR).replace(".","") + ".gpd > " + "positive_candidate_list1")

    print_run("mv refSeq_MLE_input.txt refSeq_MLE_input1.txt")

    markknownTranscripts_cmd = python_bin_foldername + "markKnownTranscripts.py " + "refSeq_MLE_input1.txt " + "positive_candidate_list1 " + "refSeq_MLE_input_marked1.txt" 
    print_run(markknownTranscripts_cmd)

    MLE_cmd = python_bin_foldername + "MLE_MT.py " + "refSeq_MLE_input_marked1.txt " + "refSeq_MLE_output1.txt " + str(Nthread) + " " + python_path
    print_run(MLE_cmd)
    reformat_cmd = python_bin_foldername + "reformat.py " + "refSeq_MLE_output1.txt > " + "refSeq_MLE_output1.txt_"
    print_run(reformat_cmd)

    maketab_cmd = python_bin_foldername + "MLEout2tab.py " + "refSeq_MLE_output1.txt_ > " + "refSeq_MLE_output1.tab"
    print_run(maketab_cmd)

    os.chdir(begin_dir)

    ####################################################

    print_run("cp " + temp_foldername + "known_LR.gpd_ref.gpd " + output_foldername + "isoform_detection.gpd")

    print_run("cp " + temp_foldername + "refSeq_MLE_output_FPR" +  str(FPR).replace(".","") + ".gpd " + output_foldername + "isoform_prediction.gpd")

    print_run("cp " + temp_foldername + "DP" +  str(FPR).replace(".","") + ".gpd " + output_foldername + "isoform.gpd")

    print_run("cp " + temp_foldername + "refSeq_MLE_output1.tab " + output_foldername + "isoform.exp")







