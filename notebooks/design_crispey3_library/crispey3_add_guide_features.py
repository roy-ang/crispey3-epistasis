#!/usr/bin/env python2
# this script adds guide features to guides table, including Azimuth guide prediction score and off-target search
# TODO: find out what's happening to some guides without XM or NM flags at all? are they not mapping?
# I suspect the is_no_insertions = (cur_num_mismatches == cur_edit_distance) check is not working for some guides in calc_guides_off_targets
import os, multiprocessing
import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from add_guide_features_functions import cal_azimuth_score, calc_guides_off_targets, add_contain_excluded_sequences_column 

def calc_guides_off_targets_unpack(arguments):
    """
    helper function to unpack arguments for calc_guides_off_targets multiprocessing
    """
    return calc_guides_off_targets(*arguments)


############################## CRISPEY LIBRARY PARAMETERS ###############################
# set library parameters
var_id_col_name = 'var_id'
guide_id_col_name = 'guide_id'
pam_seq = "GG"
excluded_seqs = ['AAAAAAAAAA',
                 'CCCCCCCCCC',
                 'GGGGGGGGGG',
                 'TTTTTTTTTT',
                 'GCATGC', # SphI cut site
                 'GGCGCGCC', # AscI cut site
                 'GCGGCCGC'] # NotI cut site
donor_length = 108
run_azimuth = True # runs Azimuth to predict guide activity -- must use python 2.7
num_of_cores = 4 # number of cores available for multiprocessing

# set guide table paths
crispey_lib_dir = os.path.expanduser("~/crispey3/library_design/Output/")
guides_table_filename = crispey_lib_dir+"all_SNPs_crispey3_GG_9bp_GUIDE.tab"
guides_with_features_table_filename = guides_table_filename.replace(".tab", "_withFeatures.tab")


# set genome fasta paths
genome_dir = os.path.expanduser("~/yeast/genomes/")
# reference genome
ref_genome_fasta_filename = genome_dir+"Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.I.fa"
# genomes to check off-targets (must have bowtie2 index)
offtarget_check_genomes_list = [("BY", ref_genome_fasta_filename), 
                           ("RM", genome_dir+"RM11-1A_SGD_2015_JRIP00000000.fsa"),
                           ("YJM", genome_dir+"YJM789_Stanford_2007_AAFW02000000_highQuality31.fsa"),
                           ("YPS", genome_dir+"YPS128.genome.fa")]

# use variant ID prefixes to determine which genomes to check off-targets
var_prefix_offtarget_check = {"ERG":[True, True, True, True],
                              "EGE":[True, True, True, False],
                              "EGD":[True, True, False, True],
                              "EGC":[True, True, False, False],
                              "EGB":[True, False, True, True],
                              "EGA":[True, False, True, False],
                              "EG9":[True, False, False, True],
                              "EG8":[True, False, False, False],
                              "GXG":[True, True, True, True],
                              "TDH":[True, False, False, False],
                              "VAL":[True, True, True, True],
                              "GXE":[True, False, False, False],
                              "HSP":[True, False, False, False],
                              "HSX":[True, False, False, False]}

# length and order of each list corresponding to variant prefix must match offtarget_check_genomes_list
# alternatively, set var_prefix_offtarget_check to None to check all guides against all offtarget_check_genomes_list


############################## LOADING WORKSPACE AND FILES ##############################
# set crispey library directory as working directory
os.chdir(crispey_lib_dir)

# load guides table
guides_table = pd.read_table(guides_table_filename, sep='\t', na_values = "")

# load reference genome
ref_genome = SeqIO.to_dict(SeqIO.parse(open(ref_genome_fasta_filename),'fasta', alphabet=generic_dna))
# set genome to uppercase
ref_genome = {chrom : sequence.upper() for chrom, sequence in ref_genome.items()}

# output dataframe
guides_with_features_df = guides_table.copy()


########################## PREDICT GUIDE EFFICACY WITH AZIMUTH ##########################
if run_azimuth:
    guides_with_features_df = cal_azimuth_score(guides_with_features_df, 
                                                    output_filename_GUIDE_withScores = "", 
                                                    guides_PAMm4p3_col_name="guide_PAM_m4p3bp")
else:
    guides_with_features_df['Azimuth'] = 1.0

    
####################### CALCULATE OFFTARGETS FOR PROVIDED GENOMES #######################
offtarget_check_args_list = []
for i in range(len(offtarget_check_genomes_list)):
    # get genome to search for off-targets
    genome_name, offtarget_check_genome_fasta_filename = offtarget_check_genomes_list[i]
       
    # load guides_df
    if var_prefix_offtarget_check is not None:
        # subset guides table by variant ID prefix if dictionary mapping provided
        var_prefix_list = [var_prefix for var_prefix, genome_check_list in var_prefix_offtarget_check.items() if genome_check_list[i]]
        guides_filter = [(v[0] in var_prefix_list) for v in guides_with_features_df[var_id_col_name].str.split('_')]
        guides_df = guides_with_features_df.loc[guides_filter,:]
    else:
        # load all guides
        guides_df = guides_with_features_df
    
    # if no guides to check off-targets in current genome, move on to next genome
    if guides_df.empty:
        continue
    
    # assemble arguments into list
    offtarget_check_args = [guides_df, offtarget_check_genome_fasta_filename, var_id_col_name, guide_id_col_name, 'guide_noPAM', pam_seq, '']
    offtarget_check_args_list.append(offtarget_check_args)

# run calc_guides_off_targets using multiprocessing
pool = multiprocessing.Pool(num_of_cores)
result = pool.map(calc_guides_off_targets_unpack, offtarget_check_args_list)

# combine offtarget results across all genomes by the variant and guide IDs in full guides table
result = [guides_with_features_df.loc[:,[var_id_col_name, guide_id_col_name]].merge(res, how='left', on=[var_id_col_name, guide_id_col_name], sort=False).fillna(0).set_index([var_id_col_name, guide_id_col_name]).astype('int') for res in result]
offtargets_df = sum(result).reset_index()
# add offtarget results to guides_with_features_df
guides_with_features_df = guides_with_features_df.merge(offtargets_df, how='left', on=[var_id_col_name, guide_id_col_name], sort=False)


################################ CHECK FOR EXCLUDED SEQS ################################
guides_with_features_df = add_contain_excluded_sequences_column(guides_with_features_df, ref_genome, excluded_seqs, donor_length)


########################## CHECK VARIANT IN NUCLEAR CHROMOSOME ##########################
guides_with_features_df['is_nuclear_chromosome'] = (guides_with_features_df['chrom'] != 'chr_mito')


##################################### WRITE TO FILE #####################################
if len(guides_with_features_table_filename) > 3:
    print("saving guides with features to: " +  guides_with_features_table_filename)
    guides_with_features_df.to_csv(guides_with_features_table_filename, sep='\t', index = False)
