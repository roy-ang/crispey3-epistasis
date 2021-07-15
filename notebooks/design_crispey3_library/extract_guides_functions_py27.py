#!/usr/bin/env python2

from __future__ import division

import sys
import os

import itertools

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import numpy as np
import math
import re

import shlex 
from subprocess import Popen, PIPE, STDOUT

import pandas as pd
import random

import vcf

import ast

import warnings
import roman

# Azimuth
#sys.path.append(os.path.expanduser("~/software/Azimuth"))
# hack - need to fix it to work with python 3

#TODO uncomment for python2_7
import azimuth.model_comparison


###############################################################################################################
# general utils
###############################################################################################################


###############################################################################################################
def is_exe(fpath):
    #print "testing:" +  fpath + ", results:" + str(os.path.isfile(fpath)) + "|" + str(os.access(fpath, os.X_OK))
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

###############################################################################################################
def which(program):
    """Replicate the UNIX which command.

    Taken verbatim from:
        stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    :program: Name of executable to test.
    :returns: Path to the program or None on failure.
    """
    fpath, program = os.path.split(program)
    if fpath:
        if is_exe(program):
            return os.path.abspath(program)
        else:
            raise ValueError(program + " is not accessible")
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return os.path.abspath(exe_file)
    
    # raise ValueError("did not find " + program + " in system path")
    return None

###############################################################################################################
# functions for guide and donor sequences design
###############################################################################################################

###############################################################################################################
def get_fasta_str(guides_df, guide_id_col_name, guide_col_name):
    """ outputs fasta file like string
    """ 
    
    fasta_str = ''.join('>' + guides_df[guide_id_col_name] + '\n' + guides_df[guide_col_name] + '\n')
    
    return(fasta_str)

##############################################################################################################
def cal_azimuth_score(guides_df, output_filename_GUIDE_withScores = "", guides_PAMm4p3_col_name="guide_PAM_m4p3bp"):
    # returns guides_df with Azimuth columns
    # output_filename_GUIDE_withScores - if not empty saves scores to file
    
    ################################################################################
    # running Azimuth score prediction (takes 20min)
    # I commented out the code that checks that the PAM sequences is NGG to allow NGA
    ################################################################################

    if (guides_df.shape[0] < 1):
        guides_df['Azimuth'] = 0
    else:
        guides_for_azimuth = np.asarray(guides_df.ix[:,[guides_PAMm4p3_col_name]]).T[0]
        CUT_POSITION = np.full_like(guides_for_azimuth, -1)
        PERCENT_PEPTIDE = np.full_like(guides_for_azimuth,-1)
        #CUT_POSITION.fill(-1)
        #PERCENT_PEPTIDE.fill(-1)
    
        #print guides_df.head()
        #print guides_for_azimuth[0:5]
        #print guides_filename
    
        #print guides_for_azimuth
        
        GUIDE_SCORES = azimuth.model_comparison.predict(guides_for_azimuth, CUT_POSITION, PERCENT_PEPTIDE)
        
        # adding the Azimuth score to the data frame
        guides_df['Azimuth'] = pd.Series(GUIDE_SCORES, index=guides_df.index)

    # write guide with scores to file
    if output_filename_GUIDE_withScores:
        guides_df.to_csv(output_filename_GUIDE_withScores, sep='\t', index = False)

    return(guides_df)

###############################################################################################################
def add_contain_excluded_sequences_column(guides_df, genome_seq, excluded_seqs, donor_length):
    """
    calculate whether a sequence around the cut site, in the guide orientation contains excluded sequences
    The input table should contain the columns:
    chrom, guide_cut_chr_pos, guide_strand 
    """
    
    seq_len_around_cut_left = int(np.floor(donor_length/2))
    seq_len_around_cut_right = int(donor_length - seq_len_around_cut_left)
    
    
    guides_df['dna_around_guide'] = ""
    
    
    print("----- start testing for excluded sequences -----")
    
    guides_df['dna_around_guide'] = ""
    
    # iterating over the table and computing the dna_around_guide
    for idx,row in guides_df.iterrows():
        
        if (idx % 20000 == 0):
            print("testing for excluded sequences in: %d" % (idx))
        
        cur_left = int(row['guide_cut_chr_pos']) - seq_len_around_cut_left
        cur_right = int(row['guide_cut_chr_pos']) + seq_len_around_cut_right
        

        cur_seq_around_cut =  genome_seq[str(row['chrom'])].seq[cur_left:cur_right]
        
        if (str(row['guide_strand']) == '-'):
            cur_seq_around_cut.reverse_complement()
        
        #row['dna_around_guide'] =  str(cur_seq_around_cut)
        #guides_df.set_value(idx,'dna_around_cut',str(cur_seq_around_cut))
        guides_df.at[idx,'dna_around_cut'] = str(cur_seq_around_cut)
    
    
    print("----- finish testing for excluded sequences -----")
    
    
    # does the DNA around the guide contains excluded sequences
    guides_df['contain_excluded_sequences'] =  guides_df['dna_around_cut'].str.contains( '|'.join(excluded_seqs) )
    
    
    return(guides_df)

###############################################################################################################
def eval_guides_off_targets (guides_df, genome_seq, guide_id_col_name, guide_col_name, pam_2nt_seq, mapping_cmd = "", input_genome_fasta_filename = ""):
    """ 
    """
    
    # parsing mapping mapping with PAM
    #guides_df['guide_map_mismatch_0'] = 0
    #guides_df['guide_map_mismatch_1'] = 0
    #guides_df['guide_map_mismatch_2'] = 0
    #guides_df['guide_map_mismatch_3'] = 0
    
    if (guides_df.shape[0] > 0):    
    
        pam_2nt_seq = Seq(pam_2nt_seq,generic_dna)
        pam_2nt_seq_revcomp = pam_2nt_seq.reverse_complement()
        
        guides_fa_str = get_fasta_str(guides_df, guide_id_col_name, guide_col_name)
        
        if len(mapping_cmd) < 1 and len(input_genome_fasta_filename) >  3:
            print("eval_guides_off_targets - composing mapping comannd using bowtie2 and the input fasta file prefix")
            bowtie_cmd = shlex.split(which('bowtie2') + ' -x ' + os.path.splitext(input_genome_fasta_filename)[0] + ' -U - -f -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 --gbar 3 --end-to-end -k 30 --no-head -t  --rdg 10,6 --rfg 10,6')
            mapping_cmd = bowtie_cmd
        
        
        print(guides_fa_str)
        
        # run mapping (bowtie2)
        #print "Start running bowtie2 for off target detection"
        mapping_pipe = Popen(mapping_cmd, stdout=PIPE, stdin=PIPE, stderr=PIPE)
        mapping_stdout = mapping_pipe.communicate(input=guides_fa_str)[0].split('\n')
        #print "Finish running bowtie2"
    
        
        
        
        #print "----- start iterating over the off target mapping -----"
        l=0
        for line in mapping_stdout:
            if (len(line)<=1):
                continue
            
            l=l+1
            if (l % 5000 == 0):
                print("Parsing mapping line: %d" % (l))
            
            #print line
            
            line_comp = line.split('\t')
            cur_guide_id  = line_comp[0]
            cur_flags     = int(line_comp[1])
            cur_chr       = line_comp[2]
            cur_left_coor = int(line_comp[3])-1 # 1- based of left most (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
            cur_read      = line_comp[9] # revcomp if mapped to the reverse strand
            
            cur_num_mismatches = None
            cur_edit_distance = None
            for comp in line_comp:
                if comp.startswith("XM:i:"):
                    cur_num_mismatches = int(comp.split(':')[2]) # XM:i:<N
                elif comp.startswith("NM:i:"):
                    cur_edit_distance = int(comp.split(':')[2]) # NM:i:<N
            
            #cur_num_mismatches = int(line_comp[14].split(':')[2]) # XM:i:<N
            #cur_edit_distance  =  int(line_comp[17].split(':')[2]) # NM:i:<N
            
            if cur_num_mismatches is None or cur_edit_distance is None:
                # guide did not return any matches, seq likely contains homopolymer, interferes with bowtie2 alignment
                # set mismatch columns to 100, thus will not pass guide filters (TODO: FIND BETTER WAY TO HANDLE SUCH GUIDES?)
                print("Warning: %s did not return a match with bowtie2. Setting guide mismatch values to 100" % cur_guide_id)
                guides_df.ix[guides_df[guide_id_col_name] == cur_guide_id, 'guide_map_mismatch_0'] = 101
                guides_df.ix[guides_df[guide_id_col_name] == cur_guide_id, 'guide_map_mismatch_1'] = 100
                guides_df.ix[guides_df[guide_id_col_name] == cur_guide_id, 'guide_map_mismatch_2'] = 100
                guides_df.ix[guides_df[guide_id_col_name] == cur_guide_id, 'guide_map_mismatch_3'] = 100
                continue
#                 raise ValueError("off target mapping: bowtie2 line does not contain XM or NM:" + line)
            
            
            ####################
            # check for pam in hit
            ####################
            is_reverse_strand = (cur_flags & 16 == 16)
            #print "id: %s, flags: %d , is reverse: %d " % (cur_guide_id, cur_flags, is_reverse_strand)
            
            if (is_reverse_strand):
                is_pam_exists = (str(genome_seq[cur_chr][(cur_left_coor-3):(cur_left_coor-1)].seq) == pam_2nt_seq_revcomp)
            else:
                is_pam_exists = (str(genome_seq[cur_chr][(cur_left_coor+21):(cur_left_coor+23)].seq) == pam_2nt_seq)
                
            
            is_no_insertions = (cur_num_mismatches == cur_edit_distance)
            
            #print ("%d , %d , %d " % (is_pam_exists, is_no_insertions, cur_num_mismatches)  )
            
            if (is_pam_exists and is_no_insertions):
                if (cur_num_mismatches == 0):
                    guides_df.ix[guides_df[guide_id_col_name] == cur_guide_id, 'guide_map_mismatch_0'] += 1
                elif (cur_num_mismatches == 1):
                    guides_df.ix[guides_df[guide_id_col_name] == cur_guide_id, 'guide_map_mismatch_1'] += 1
                elif (cur_num_mismatches == 2):
                    guides_df.ix[guides_df[guide_id_col_name] == cur_guide_id, 'guide_map_mismatch_2'] += 1
                elif (cur_num_mismatches == 3):
                    guides_df.ix[guides_df[guide_id_col_name] == cur_guide_id, 'guide_map_mismatch_3'] += 1
     
     
        #print "----- finish iterating over the off target mapping -----"
        
        # removing mapping on target
        guides_df['guide_map_mismatch_0'] =  guides_df['guide_map_mismatch_0'] - 1
        
    return(guides_df)


###############################################################################################################
# high level function for guide and donor sequence design
###############################################################################################################

###############################################################################################################
def cal_guide_features_multiGenomeOffTargetCheck(input_guides_table_filename, input_genome_fasta_filename,
                                                 offtarget_check_genomes_fastas, output_guides_with_features_table_filename = "", 
                                                 PAM_seq = 'GA', excluded_seqs = ['TTTTT'], donor_length = 100,
                                                 BOWTIE_exe = which('bowtie2'), run_azimuth = False):
    """
    NEEDS TO BE FIXED:
    set up function(s) to determine which variants need to be checked with which genomes.
    a higher-level function (crispr_add_guide_features) can first subset the guides table according to which variants need
    to be checked against which genomes. For each genome, cal_guide_features can be run with a specific subset guides table.
    Once all genomes are checked, table can be recombined and off-target counts added together.
    Supply a dictionary of variant ID prefixes as keys, that each correspond to a list of booleans, length equal to the number 
    of genomes to check for offtargets. This dictionary will be first parsed to determine the subsetting of guide table for 
    each tested genome
    
    or just ditch the higher level script and write a single python script with variables stored within in... no need to get fancy
    
    
    
    Calculate guides features

    input_genome_fasta_filename - must have a bowtie2 index build 
    input_guides_table_filename - must contain: guide_id, guide_noPAM, guide_PAM_m4p3bp
    
    Adds to the table
    
    1. Azimuth

    Off target:
    2. guide_map_mismatch_0 
    3. guide_map_mismatch_1 
    4. guide_map_mismatch_2 
    5. guide_map_mismatch_3

    6. contain_excluded_sequences
    7. is_nuclear_chromosome
    8. dna_around_guide
    
    Off target analysis can be conducted for multiple genomes, reports a number to guide_map_mismatch for each genome analyzed

    """

    
    # fixing the mapping command to be the default bowtie
    if not BOWTIE_exe or len(BOWTIE_exe) < 3:
        
        print("input bowtie2 is null, using which function to find path")
        BOWTIE_exe = which('bowtie2')
        
        if not BOWTIE_exe:
            raise ValueError("bowtie2 is unknown (which returns None, make sure unix >>which bowtie2)")
        print('###')
        


    # defining the bowtie mapping cmd 
    # (this has been disabled to allow mapping cmd to be set for each genome used for offtarget checks)
    #bowtie_cmd = shlex.split(BOWTIE_exe + ' -x ' + os.path.splitext(input_genome_fasta_filename)[0] + ' -U - -f -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 --gbar 3 --end-to-end -k 30 --no-head -t  --rdg 10,6 --rfg 10,6')
    #mapping_cmd = bowtie_cmd
    
    # TODO remove this printing 
    #print("Using this command for mapping off targets: " + BOWTIE_exe + ' -x ' + os.path.splitext(input_genome_fasta_filename)[0] + ' -U - -f -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 --gbar 3 --end-to-end -k 30 --no-head -t  --rdg 10,6 --rfg 10,6')

    # leave mapping cmd input empty, allow eval_guides_off_targets to assemble mapping cmd using genome_filename
    mapping_cmd = ""
    
    # load guide df
    guides_df = pd.read_table(input_guides_table_filename, sep='\t', na_values = "")

    # loading genome fasta file
    genome_seq = SeqIO.to_dict(SeqIO.parse(open(input_genome_fasta_filename),'fasta', alphabet=generic_dna))
 
    # copying the original df matrix TODO removed the copying since uplaoded from file
    out_guides_with_features_df = guides_df #.copy()
    
    # calculating Azimuth score
    if run_azimuth:
        out_guides_with_features_df = cal_azimuth_score(out_guides_with_features_df, 
                                                        output_filename_GUIDE_withScores = "", 
                                                        guides_PAMm4p3_col_name="guide_PAM_m4p3bp")
    else:
        out_guides_with_features_df['Azimuth'] = 1.0
    
    #DEBUG
    #print "XXXX"
    
    # calculate off targets
    out_guides_with_features_df['guide_map_mismatch_0'] = 0
    out_guides_with_features_df['guide_map_mismatch_1'] = 0
    out_guides_with_features_df['guide_map_mismatch_2'] = 0
    out_guides_with_features_df['guide_map_mismatch_3'] = 0
    for genome_filename in offtarget_check_genomes_fastas:
        print("Check off targets using genome", genome_filename)
        genome = SeqIO.to_dict(SeqIO.parse(open(genome_filename),'fasta', alphabet=generic_dna))
        out_guides_with_features_df = eval_guides_off_targets(out_guides_with_features_df, genome, 
                                                              guide_id_col_name = 'guide_id', guide_col_name = 'guide_noPAM', 
                                                              pam_2nt_seq = PAM_seq, mapping_cmd = mapping_cmd,
                                                              input_genome_fasta_filename = genome_filename)
    
    #DEBUG
    #print "after eval_guides_off_targets"
    #print out_guides_with_features_df
    
    # calculating if containing excluded sequences    
    out_guides_with_features_df = add_contain_excluded_sequences_column(out_guides_with_features_df, genome_seq, excluded_seqs, donor_length)

    # calculate if in the major nuclear chromosome
    out_guides_with_features_df['is_nuclear_chromosome'] = (out_guides_with_features_df['chrom'] != 'chr_mito')
    
    if len(output_guides_with_features_table_filename) > 3:
            print("saving guides with features to: " +  output_guides_with_features_table_filename)
            out_guides_with_features_df.to_csv(output_guides_with_features_table_filename, sep='\t', index = False)
    
    # returning the updated data frame
    return(out_guides_with_features_df)

