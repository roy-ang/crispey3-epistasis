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

# from gff2dataframe import sgd_gff2dataframe

# Azimuth
#sys.path.append(os.path.expanduser("~/software/Azimuth"))
# hack - need to fix it to work with python 3

#TODO uncomment for python2_7
#import azimuth.model_comparison


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
def file_len(fname):
    # not including lines that start with #
    i=0
    with open(fname) as f:
        comment_char = '#'
        comment_lines = 0
        for i, l in enumerate(f,1):
            if l[0] == comment_char:
                comment_lines = comment_lines + 1
        print("# comments lines: %d" % (comment_lines))
        i = i - comment_lines
    return i



###############################################################################################################
# roman
###############################################################################################################

def int2roman(num):

    roman = OrderedDict()
    roman[1000] = "M"
    roman[900] = "CM"
    roman[500] = "D"
    roman[400] = "CD"
    roman[100] = "C"
    roman[90] = "XC"
    roman[50] = "L"
    roman[40] = "XL"
    roman[10] = "X"
    roman[9] = "IX"
    roman[5] = "V"
    roman[4] = "IV"
    roman[1] = "I"

    def roman_num(num):
        for r in roman.keys():
            x, y = divmod(num, r)
            yield roman[r] * x
            num -= (r * x)
            if num > 0:
                roman_num(num)
            else:
                break
              
    return("".join([d for d in roman_num(num)]))


def yeast_chr_numeric2roman(chr_numeric_str):
    return('chr' + int2roman(int(chr_numeric_str.replace('chr',''))))

###############################################################################################################
# VEP annotations
###############################################################################################################

def get_genes_len(input_gff_filename):
    with open(input_gff_filename) as f:
        gff_lines = f.readlines()
    gff_lines = [l for l in gff_lines if len(l.split('\t'))>=5 ]    
    gff_lines = [l for l in gff_lines if l.split('\t')[2] == 'gene' ]
    genes_ids = [l.split('ID=')[1].split(';Name=')[0] for l in gff_lines]
    genes_lens = [ int(l.split('\t')[4]) - int(l.split('\t')[3]) + 1 for l in gff_lines]

    gene_le_df = pd.DataFrame({'Gene' : genes_ids, 'len_bp' : genes_lens})
    gene_le_df.set_index(keys='Gene', drop=False, inplace=True)
    gene_le_df.index.name=''
    gene_le_df.drop_duplicates(inplace=True)

    return(gene_le_df)


###############################################################################################################

def get_cdss_len(input_gff_filename):
    with open(input_gff_filename) as f:
        gff_lines = f.readlines()
    gff_lines = [l for l in gff_lines if len(l.split('\t'))>=5 ]
    gff_lines = [l for l in gff_lines if l.split('\t')[2] == 'CDS' ]
    genes_ids = [l.split('Parent=')[1].split(';Name=')[0] for l in gff_lines]
    genes_lens = [ int(l.split('\t')[4]) - int(l.split('\t')[3]) + 1 for l in gff_lines]

    gene_le_df = pd.DataFrame({'Gene' : genes_ids, 'len_bp' : genes_lens})
    gene_le_df.set_index(keys='Gene', drop=False, inplace=True)
    gene_le_df.index.name=''
    gene_le_df = gene_le_df.groupby('Gene').sum()
    gene_le_df['Gene'] = gene_le_df.index
    gene_le_df.index.name=''
    gene_le_df.drop_duplicates(inplace=True)

    return(gene_le_df)

###############################################################################################################

def get_dubious_genes(input_gff_filename):
    with open(input_gff_filename) as f:
        gff_lines = f.readlines()
    gff_lines = [l for l in gff_lines if len(l.split('\t'))>=5 ]
    gff_lines = [l for l in gff_lines if (l.split('\t')[2] == 'gene') & (len(re.findall('dubious', l, flags=re.I))>0)]
    genes_ids = [l.split('ID=')[1].split(';Name=')[0] for l in gff_lines]
    
    return genes_ids


###############################################################################################################

def annotate_variants_by_VEPoutput(vep_input_vcf_filename, vep_output_filename, input_gff_filename, annotated_res_output_filename, id_colname = 'var_id'):
    cdss_len_df = get_cdss_len(input_gff_filename)
    genes_len_df = get_genes_len(input_gff_filename)
    dubious_genes = get_dubious_genes(input_gff_filename)
    
    vep_df = pd.read_table(vep_output_filename, sep='\t', na_values = "-", low_memory=False)
    vep_df = vep_df.loc[~vep_df.Gene.isin(dubious_genes),]
    vep_df.rename(columns={"#Uploaded_variation" : id_colname}, inplace=True)
    
    # prepare annotated variants df using VCF supplied to VEP
    variants_annotated_df = pd.Series([rec for rec in vcf.Reader(filename=vep_input_vcf_filename)])
    var_id_list = []
    assemble = []
    columns = ['CHROM', 'POS', 'REF', 'ALT', 'AC', 'AN' ,'AF']
    for rec in variants_annotated_df:
        var_id_list.append(rec.ID)
        assemble.append((rec.CHROM, rec.POS, rec.REF, rec.ALT[0], rec.INFO['AC'][0], rec.INFO['AN'], rec.INFO['AF'][0]))
    variants_annotated_df = pd.DataFrame.from_records(assemble, index=var_id_list, columns=columns)
    variants_annotated_df.index.name = id_colname
    
    # parse VEP output text file by ID for annotations
    entry_cnt = 0
    for cur_id, anno in vep_df.groupby(id_colname):
        entry_cnt +=1
        if (entry_cnt % 5000 == 0):
            print("parsing entry: %d"% (entry_cnt))
        
        # all up/down stream genes (VEP default search range is 5kb from variant position)
        cur_upstream_all_str =  '|'.join([str(int(d)) for d in anno.DISTANCE[anno.Consequence.isin(['upstream_gene_variant' ])].values if not np.isnan(d)])
        cur_downstream_all_str =  '|'.join([str(int(d)) for d in anno.DISTANCE[anno.Consequence.isin(['downstream_gene_variant' ])].values if not np.isnan(d)])
        
        variants_annotated_df.loc[cur_id, 'upstream_distance_str'] = cur_upstream_all_str
        variants_annotated_df.loc[cur_id, 'downstream_distance_str'] = cur_downstream_all_str

        is_up_down_stream = anno.Consequence.isin(['upstream_gene_variant' ,'downstream_gene_variant'])
        
        # identify intergenic variants
        if np.all(is_up_down_stream):
            left_genes = ( ((anno.Consequence =='upstream_gene_variant') & (anno.STRAND<0)) |
                          ((anno.Consequence =='downstream_gene_variant') & (anno.STRAND>0)) )

            right_genes = ( ((anno.Consequence =='upstream_gene_variant') & (anno.STRAND>0)) |
                              ((anno.Consequence =='downstream_gene_variant') & (anno.STRAND<0)) )

            if (left_genes.sum() > 0):
                ind_of_closest_left = anno.loc[left_genes,'DISTANCE'].idxmin()
                
                variants_annotated_df.loc[cur_id, 'closest_gene1_Gene_Name'] = anno['SYMBOL'][ind_of_closest_left]
                variants_annotated_df.loc[cur_id, 'closest_gene1_Gene_ID'] = anno['Gene'][ind_of_closest_left]
                variants_annotated_df.loc[cur_id, 'closest_gene1_Annotation'] = anno['Consequence'][ind_of_closest_left]
                variants_annotated_df.loc[cur_id, 'closest_gene1_Distance'] = anno['DISTANCE'][ind_of_closest_left]

            if (right_genes.sum() > 0):
                ind_of_closest_right = anno.loc[right_genes,'DISTANCE'].idxmin()
                
                variants_annotated_df.loc[cur_id, 'closest_gene2_Gene_Name'] = anno['SYMBOL'][ind_of_closest_right]
                variants_annotated_df.loc[cur_id, 'closest_gene2_Gene_ID'] = anno['Gene'][ind_of_closest_right]
                variants_annotated_df.loc[cur_id, 'closest_gene2_Annotation'] = anno['Consequence'][ind_of_closest_right]
                variants_annotated_df.loc[cur_id, 'closest_gene2_Distance'] = anno['DISTANCE'][ind_of_closest_right]

            # determine type of noncoding variant
            # Caution: variants without genes on one side will be labelled as "intergenic"
            cur_region = 'intergenic'
            if ((left_genes.sum() > 0) and (right_genes.sum() > 0)): 
                if (anno['Consequence'][ind_of_closest_left] == 'upstream_gene_variant'):
                    if (anno['Consequence'][ind_of_closest_right] == 'upstream_gene_variant'):
                        cur_region = 'bidirectional_promoter'
                    else:
                        cur_region = 'unidirectional_promoter'
                else:
                    if (anno['Consequence'][ind_of_closest_right] == 'upstream_gene_variant'):
                        cur_region = 'unidirectional_promoter'

            variants_annotated_df.loc[cur_id, 'region'] = cur_region

            # index of closest gene
            sel_row = anno.DISTANCE.idxmin()
        
        # for variants that land in a gene, select annotation with highest impact
        else:
            anno = anno.loc[~is_up_down_stream,:]
            if (anno.IMPACT == 'HIGH').sum()>0:
                anno = anno.loc[anno.IMPACT == 'HIGH',:]
            elif (anno.IMPACT == 'MODERATE').sum()>0:
                anno = anno.loc[anno.IMPACT == 'MODERATE',:]
            elif (anno.IMPACT == 'LOW').sum()>0:
                anno = anno.loc[anno.IMPACT == 'LOW',:]
            elif (anno.IMPACT == 'MODIFIER').sum()>0:
                anno = anno.loc[anno.IMPACT == 'MODIFIER',:] #NOTE: includes introns, noncoding transcripts within length of gene
            
            else:
                raise ValueError(anno)

            sel_row = anno.index[0]
            variants_annotated_df.loc[cur_id, 'region'] = anno['Consequence'][sel_row]

        # closest annotation recorded
        variants_annotated_df.loc[cur_id, 'Annotation'] = anno['Consequence'][sel_row]
        variants_annotated_df.loc[cur_id, 'Annotation_Impact'] = anno['IMPACT'][sel_row]
        variants_annotated_df.loc[cur_id, 'Gene_Name'] = anno['SYMBOL'][sel_row]
        variants_annotated_df.loc[cur_id, 'Gene_ID'] = anno['Gene'][sel_row]
        variants_annotated_df.loc[cur_id, 'Transcript_BioType'] = anno['BIOTYPE'][sel_row]
        variants_annotated_df.loc[cur_id, 'Existing_variation'] = anno['Existing_variation'][sel_row]
        variants_annotated_df.loc[cur_id, 'HGVSc'] = anno['HGVSc'][sel_row]
        variants_annotated_df.loc[cur_id, 'HGVSp'] = anno['HGVSp'][sel_row]
        variants_annotated_df.loc[cur_id, 'SWISSPROT'] = anno['SWISSPROT'][sel_row]
        variants_annotated_df.loc[cur_id, 'UNIPARC'] = anno['UNIPARC'][sel_row]
        variants_annotated_df.loc[cur_id, 'BLOSUM62'] = anno['BLOSUM62'][sel_row]
        variants_annotated_df.loc[cur_id, 'HGVSp'] = anno['HGVSp'][sel_row]
#         variants_annotated_df.set_value(index=cur_id,col='VEP_CSN',value=anno['CSN'][sel_row])
        
        # store cDNA and CDS positions
        cdna_pos = anno['cDNA_position'][sel_row]
        cds_pos = anno['CDS_position'][sel_row]
    
        # add additional info for missense and synonymous variants
        if (anno['Consequence'][sel_row] in ['synonymous_variant','missense_variant']):
            # because variants were extracted from 1011 genomes GVCF, some SNPs have trailing bases.
            # Calculate the correct cDNA and CDS positions if so 
            if '-' in cdna_pos:
                cdna_pos = cdna_pos.split('-')[0]
            if '-' in cds_pos:
                cds_pos = cds_pos.split('-')[0]
                
            
            if anno['Gene'][sel_row] in genes_len_df.index:
                cur_gene_len = int(genes_len_df.loc[anno['Gene'][sel_row],'len_bp'])
                variants_annotated_df.loc[cur_id, 'cDNA_len'] = cur_gene_len
                variants_annotated_df.loc[cur_id, 'cDNA_frac'] = int(cdna_pos)/cur_gene_len
            
            if anno['Gene'][sel_row] in cdss_len_df.index:
                cur_cds_len = int(cdss_len_df.loc[anno['Gene'][sel_row],'len_bp'])
                variants_annotated_df.loc[cur_id, 'CDS_len'] = cur_cds_len
                variants_annotated_df.loc[cur_id, 'CDS_frac'] = int(cds_pos)/cur_cds_len
                variants_annotated_df.loc[cur_id, 'AA_len'] = cur_cds_len/3
                variants_annotated_df.loc[cur_id, 'AA_frac'] = int(cds_pos)/cur_cds_len
                
            variants_annotated_df.loc[cur_id, 'cDNA_pos'] = cdna_pos
            variants_annotated_df.loc[cur_id, 'CDS_pos'] = cds_pos
            variants_annotated_df.loc[cur_id, 'AA_pos'] = anno['Protein_position'][sel_row]
            variants_annotated_df.loc[cur_id, 'Amino_acids'] = anno['Amino_acids'][sel_row]
            variants_annotated_df.loc[cur_id, 'Codons'] = anno['Codons'][sel_row]
            variants_annotated_df.loc[cur_id, 'DOMAINS'] = anno['DOMAINS'][sel_row]

    variants_annotated_df.to_csv(annotated_res_output_filename, sep='\t', header=True, index=True)
    
    return(variants_annotated_df)


###############################################################################################################

def parse_VEF_output_fromSNP_old2018(vep_input_df, vep_output_filename, input_gff_filename,
                    annotated_res_output_filename, id_colname = 'oligo_seq'):
    
    
    cdss_len_df = get_cdss_len(input_gff_filename)
    genes_len_df = get_genes_len(input_gff_filename)
    
    vep_df = pd.read_table(vep_output_filename, sep='\t', na_values = "-")
    dubious_genes = get_dubious_genes(input_gff_filename)
    vep_df = vep_df.loc[~vep_df.Feature.isin(dubious_genes),]
    vep_df.rename(columns={"#Uploaded_variation" : id_colname}, inplace=True)
    
    annotate_df = create_annot_for_VEP_fromSNP(vep_input_df, id_colname)
    
    vep_df_group_by_id = vep_df.groupby(id_colname)
    
    entry_cnt = 0
    for cur_id, oligo_anns in vep_df_group_by_id:
        entry_cnt +=1
        
        if (entry_cnt % 5000 == 0):
            print("parsing entry: %d"% (entry_cnt))
        

        # all up/down stream genes up to 5k
        cur_upstream_all_str =  '|'.join([str(int(d)) for d in oligo_anns.DISTANCE[oligo_anns.Consequence.isin(['upstream_gene_variant' ])].values if not np.isnan(d)])
        cur_downstream_all_str =  '|'.join([str(int(d)) for d in oligo_anns.DISTANCE[oligo_anns.Consequence.isin(['downstream_gene_variant' ])].values if not np.isnan(d)])

        annotate_df.set_value(index=cur_id,col='SNPEff_upstream_distance_str',value=cur_upstream_all_str)
        annotate_df.set_value(index=cur_id,col='SNPEff_downstream_distance_str',value=cur_downstream_all_str)

        is_up_down_stream_tf = oligo_anns.Consequence.isin(['upstream_gene_variant' ,'downstream_gene_variant'])

        if np.all(is_up_down_stream_tf):
            # intergenic
            left_genes_tf = ( ((oligo_anns.Consequence =='upstream_gene_variant') & (oligo_anns.STRAND<0)) |
                          ((oligo_anns.Consequence =='downstream_gene_variant') & (oligo_anns.STRAND>0)) )

            right_genes_tf = ( ((oligo_anns.Consequence =='upstream_gene_variant') & (oligo_anns.STRAND>0)) |
                              ((oligo_anns.Consequence =='downstream_gene_variant') & (oligo_anns.STRAND<0)) )

            if (left_genes_tf.sum() > 0):
                ind_of_closest_left = np.argmin(oligo_anns.loc[left_genes_tf,'DISTANCE'])

                annotate_df.set_value(index=cur_id,col='SNPEff_closets_gene1_Gene_Name',value=oligo_anns['SYMBOL'][ind_of_closest_left])
                annotate_df.set_value(index=cur_id,col='SNPEff_closets_gene1_Gene_ID',value=oligo_anns['Gene'][ind_of_closest_left])
                annotate_df.set_value(index=cur_id,col='SNPEff_closets_gene1_Annotation',value=oligo_anns['Consequence'][ind_of_closest_left])
                annotate_df.set_value(index=cur_id,col='SNPEff_closets_gene1_Distance',value=oligo_anns['DISTANCE'][ind_of_closest_left])

            if (right_genes_tf.sum() > 0):
                ind_of_closest_right = np.argmin(oligo_anns.loc[right_genes_tf,'DISTANCE'])

                annotate_df.set_value(index=cur_id,col='SNPEff_closets_gene2_Gene_Name',value=oligo_anns['SYMBOL'][ind_of_closest_right])
                annotate_df.set_value(index=cur_id,col='SNPEff_closets_gene2_Gene_ID',value=oligo_anns['Gene'][ind_of_closest_right])
                annotate_df.set_value(index=cur_id,col='SNPEff_closets_gene2_Annotation',value=oligo_anns['Consequence'][ind_of_closest_right])
                annotate_df.set_value(index=cur_id,col='SNPEff_closets_gene2_Distance',value=oligo_anns['DISTANCE'][ind_of_closest_right])


            # type of intergenic region
            cur_region = 'intergenic'
            if ((left_genes_tf.sum() > 0) and (right_genes_tf.sum() > 0)):
                if (oligo_anns['Consequence'][ind_of_closest_left] == 'upstream_gene_variant'):
                    if (oligo_anns['Consequence'][ind_of_closest_right] == 'upstream_gene_variant'):
                        cur_region = 'divergent promoter'
                    else:
                        cur_region = 'unidirectional promoter'
                else:
                    if (oligo_anns['Consequence'][ind_of_closest_right] == 'upstream_gene_variant'):
                        cur_region = 'unidirectional promoter'


            annotate_df.set_value(index=cur_id,col='SNPEff_region',value=cur_region)

            # index of closest gene
            sel_row = np.argmin(oligo_anns.DISTANCE)
        else:
            oligo_anns = oligo_anns.loc[~is_up_down_stream_tf,:]
            if (oligo_anns.IMPACT == 'HIGH').sum()>0:
                oligo_anns = oligo_anns.loc[oligo_anns.IMPACT == 'HIGH',:]
            elif (oligo_anns.IMPACT == 'MODERATE').sum()>0:
                oligo_anns = oligo_anns.loc[oligo_anns.IMPACT == 'MODERATE',:]
            elif (oligo_anns.IMPACT == 'MODIFIER').sum()>0:
                oligo_anns = oligo_anns.loc[oligo_anns.IMPACT == 'MODIFIER',:]
            elif (oligo_anns.IMPACT == 'LOW').sum()>0:
                oligo_anns = oligo_anns.loc[oligo_anns.IMPACT == 'LOW',:]
            else:
                raise ValueError(oligo_anns)

            sel_row = oligo_anns.index[0]
            annotate_df.set_value(index=cur_id,col='SNPEff_region',value=oligo_anns['Consequence'][sel_row])

        # closest annotation
        annotate_df.set_value(index=cur_id,col='SNPEff_Annotation',value=oligo_anns['Consequence'][sel_row])
        annotate_df.set_value(index=cur_id,col='SNPEff_Annotation_Impact',value=oligo_anns['IMPACT'][sel_row])
        annotate_df.set_value(index=cur_id,col='SNPEff_Gene_Name',value=oligo_anns['SYMBOL'][sel_row])
        annotate_df.set_value(index=cur_id,col='SNPEff_Gene_ID',value=oligo_anns['Gene'][sel_row])
        annotate_df.set_value(index=cur_id,col='SNPEff_Transcript_BioType',value=oligo_anns['BIOTYPE'][sel_row])
        annotate_df.set_value(index=cur_id,col='VEP_Existing_variation',value=oligo_anns['Existing_variation'][sel_row])
        annotate_df.set_value(index=cur_id,col='VEP_HGVSc',value=oligo_anns['HGVSc'][sel_row])
        annotate_df.set_value(index=cur_id,col='VEP_HGVSp',value=oligo_anns['HGVSp'][sel_row])

        annotate_df.set_value(index=cur_id,col='VEP_SWISSPROT',value=oligo_anns['SWISSPROT'][sel_row])
        annotate_df.set_value(index=cur_id,col='VEP_UNIPARC',value=oligo_anns['UNIPARC'][sel_row])
        
        annotate_df.set_value(index=cur_id,col='VEP_BLOSUM62',value=oligo_anns['BLOSUM62'][sel_row])
        annotate_df.set_value(index=cur_id,col='VEP_CSN',value=oligo_anns['CSN'][sel_row])

    

        if (oligo_anns['Consequence'][sel_row] in ['synonymous_variant','missense_variant']):

            if oligo_anns['Gene'][sel_row] in genes_len_df.index:
                #print('------')
                #print(genes_len_df.loc[oligo_anns['Gene'][sel_row],'len_bp'])
                #print('------')
                cur_gene_len = int(genes_len_df.loc[oligo_anns['Gene'][sel_row],'len_bp'])
                annotate_df.set_value(index=cur_id,col='SNPEff_cDNA_len',value=cur_gene_len)
                annotate_df.set_value(index=cur_id,col='SNPEff_cDNA_frac',value=int(oligo_anns['cDNA_position'][sel_row])/cur_gene_len)
                
            
            if oligo_anns['Gene'][sel_row] in cdss_len_df.index:
                #print('------')
                #print(cdss_len_df.loc[oligo_anns['Gene'][sel_row],'len_bp'])
                #print('------')
                cur_cds_len = int(cdss_len_df.loc[oligo_anns['Gene'][sel_row],'len_bp'])
                annotate_df.set_value(index=cur_id,col='SNPEff_CDS_len',value=cur_cds_len)
                annotate_df.set_value(index=cur_id,col='SNPEff_CDS_frac',value=int(oligo_anns['CDS_position'][sel_row])/cur_cds_len)
                
                annotate_df.set_value(index=cur_id,col='SNPEff_AA_len',value=cur_cds_len/3)
                annotate_df.set_value(index=cur_id,col='SNPEff_AA_frac',value=int(oligo_anns['CDS_position'][sel_row])/cur_cds_len)

         
            annotate_df.set_value(index=cur_id,col='SNPEff_cDNA_pos',value=oligo_anns['cDNA_position'][sel_row])
            annotate_df.set_value(index=cur_id,col='SNPEff_CDS_pos',value=oligo_anns['CDS_position'][sel_row])
            annotate_df.set_value(index=cur_id,col='SNPEff_AA_pos',value=oligo_anns['Protein_position'][sel_row])

            annotate_df.set_value(index=cur_id,col='VEP_Amino_acids',value=oligo_anns['Amino_acids'][sel_row])
            annotate_df.set_value(index=cur_id,col='VEP_Codons',value=oligo_anns['Codons'][sel_row])
            annotate_df.set_value(index=cur_id,col='VEP_DOMAINS',value=oligo_anns['DOMAINS'][sel_row])


    annotate_df.to_csv(annotated_res_output_filename, sep='\t', index=False)
    
    return(annotate_df)



def prepare_VEP_input_from_SNPs_df_old2018(snps_df_filename, vep_input_filename):
    snps_df = pd.read_table(snps_df_filename, sep='\t', na_values = "")
    annotate_df = snps_df.loc[snps_df.chrom != 'chr_mito', ['oligo_seq', 'var_id', 'chrom', 'SNP_chr_pos', 'REF', 'ALT', 'donor_info_str'] ].copy()

    annotate_df['ALT'] = annotate_df['ALT'].str.replace('[','').str.replace(']','')
    annotate_df['coor_start'] = annotate_df['SNP_chr_pos']
    annotate_df['coor_end'] = annotate_df['SNP_chr_pos'] + annotate_df['REF'].apply(lambda s: len(s)) - 1
    if annotate_df.chrom.str.isalpha().all():
        annotate_df['chrom_roman'] = annotate_df.chrom.str.replace('chr','')
    else:
        annotate_df['chrom_roman'] = annotate_df.chrom.apply(yeast_chr_numeric2roman).str.replace('chr','')
    annotate_df['strand'] = '+'

    annotate_df.loc[annotate_df.REF == '','REF'] = '-'
    annotate_df.loc[annotate_df.ALT == '','ALT'] = '-'
    annotate_df['allele_edit'] = annotate_df['REF'] + '/' + annotate_df['ALT']


    vep_input_df = annotate_df.loc[:,
                                   ['chrom_roman','coor_start', 'coor_end',
                                    'allele_edit', 'strand', 'var_id']]

    vep_input_df = vep_input_df.loc[np.isfinite(vep_input_df.coor_end),:]
    vep_input_df.ix[:,'coor_start'] = vep_input_df.coor_start.astype(int)
    vep_input_df.ix[:,'coor_end'] = vep_input_df.coor_end.astype(int)
    
    # writing VCF file
    vep_input_df.to_csv(vep_input_filename,sep='\t',index=False, header=False)

    return(vep_input_df)


###############################################################################################################
# general sequence utils
###############################################################################################################


###############################################################################################################
def get_random_dna_nts(length):
    return( Seq(''.join([random.choice("ACGT") for i in xrange(length)]), generic_dna) )

###############################################################################################################
def get_all_random_dna_nts(length):
    dna_bases=['A','C','G','T']
    return([''.join(p) for p in itertools.product(dna_bases, repeat=length)])

###############################################################################################################
def get_K_random_dna_nts(length,K):
    dna_bases=['A','C','G','T']
    
    # if the length is very long - create each sequnce separately and make it is different
    # else produce all sequences and select from them
    num_all_seqs = 4**length
    if (K * 100 < num_all_seqs):
        i = 0
        all_seqs = []
        # the i is a safty garde for infinate loop
        while ( (i < 100 * K) and (len(all_seqs) < K) ):
            i += 1
            cur_seq = get_random_dna_nts(length)
            if not (cur_seq in all_seqs):
                all_seqs.append(cur_seq)
                
    else:
        all_seqs = [''.join(p) for p in itertools.product(dna_bases, repeat=length)]
    
    
    if (K>len(all_seqs)):
            warnings.warn("get_K_random_dna_nts: number of optional sequences %d is smaller than the requested %d" % (len(all_seqs), K) )
       
        
    if len(all_seqs) == 0:
        return None
    else:
        return( (random.sample(all_seqs,min(K,len(all_seqs))), K<=len(all_seqs) ) )
       

###############################################################################################################
# functions for guide and donor sequences design
###############################################################################################################

##############################################################################################################
def write_var_id(record, var_id_prefix, i):
    """
    takes VCF record, var_id_prefix string and number i as input
    returns a string of var_id to be printed
    """
    if var_id_prefix == '':
        return record.ID
    else:
        return var_id_prefix + str(i)


###############################################################################################################
def get_fasta_str(guides_df, guide_id_col_name, guide_col_name):
    """ outputs fasta file like string
    """ 
    
    fasta_str = ''.join('>' + guides_df[guide_id_col_name] + '\n' + guides_df[guide_col_name] + '\n')
    
    return(fasta_str)

###############################################################################################################
def get_min_distance_to_pam(seq_record, pam_seqs):
    
    distance_to_pam = float('nan')
    minimal_dist_pam = ""
    for pam in pam_seqs:    
        cur_dist = seq_record.seq.find(pam)
        if (cur_dist > -1 and (math.isnan(distance_to_pam) or abs(cur_dist) < abs(distance_to_pam)) ):
            distance_to_pam = cur_dist
            minimal_dist_pam = pam
        #print "PAM  search %s (%d)" % (pam,cur_dist)
    
    return(distance_to_pam, minimal_dist_pam)
    

###############################################################################################################
def get_all_guide_SNP_position(seq, pam, alt):
    # numbering is   XXX  NG[G|A]
    #              -3-2-1 012
    # SNP in zero is removed
    
    # returns a list of all the positions of the SNP in the PAM
    cur_pam_dists = np.array([m.start() for m in re.finditer("(?="+pam+")", str(seq))])
    
    #print(str(seq))
    
    #print "XXX"
    #print cur_pam_dists
    #print "XXX"
    
    # positions in the guide (-k:-1)
    cur_pam_dists[cur_pam_dists > 2] = 2 - cur_pam_dists[cur_pam_dists > 2]
    
    # in case SNP is in the N of the PAM (NG[G|A]) - delte
    cur_pam_dists = cur_pam_dists[cur_pam_dists!=2]
    
    # in case SNP is in the PAM third position (2) make sure the alternative is not a valid PAM
    #if (len(alt)>1 or (alt[0] != "G" and alt[0] != "A")):
    #    cur_pam_dists[cur_pam_dists==0] = 2
    #else:
    #    cur_pam_dists = cur_pam_dists[cur_pam_dists!=0]
    
    # do not check if alternative is A/G
    cur_pam_dists[cur_pam_dists==0] = 2
    
    cur_pam_dists[cur_pam_dists==1] = 1
    
    #print "YYY"
    #print cur_pam_dists
    #print "YYY"
    
    return(cur_pam_dists.astype(int))
    

###############################################################################################################
def get_best_guide_SNP_position(seq, pam_seqs, alt = "N"):
    # seq is a string starting at -1 the SNP position
    # in no guide is found returns NaN 
    # TODO DEBUG hard coding for PAM GG / GA checking the SNP

    best_guide_SNP_position = np.nan
    best_pam = ""
    
    all_pam_SNP_positions = np.array([])
    all_pam_per_position = []
    
    for pam in pam_seqs:
        cur_pam_positions = get_all_guide_SNP_position(seq, pam, alt)
        # concatenating the results for different PAMs
        all_pam_SNP_positions = np.concatenate((all_pam_SNP_positions, cur_pam_positions),axis=0)
        all_pam_per_position = all_pam_per_position + [pam]*len(cur_pam_positions)
        #print "ALL:" + str(all_pam_SNP_positions)
        #print "ALL pam:" + str(all_pam_per_position)
    
    
    # selecting the the best guide
    
    if (len(all_pam_SNP_positions)>0):
        
        # if there are SNPs in the PAM
        if (any(all_pam_SNP_positions>0)):
            pos_inds = np.array([i for i in range(0,len(all_pam_SNP_positions)) if all_pam_SNP_positions[i]>0])
            best_ind = pos_inds[all_pam_SNP_positions[pos_inds]==all_pam_SNP_positions[pos_inds].min()][0]
            #print "Best index in PAM: %d" % (best_ind)
        else:
            best_ind = all_pam_SNP_positions.argmax()
            #print "Best index in guide: %d" % (best_ind)
        
        
        best_guide_SNP_position =   all_pam_SNP_positions[best_ind]
        best_pam = all_pam_per_position[best_ind]
        
    #else:
    #    print "No guide found"
    

     
        
#        cur_dist = seq_record.seq.find(pam)
#        if (cur_dist > -1 and (math.isnan(distance_to_pam) or abs(cur_dist) < abs(distance_to_pam)) ):
#            distance_to_pam = cur_dist
#            minimal_dist_pam = pam
#        #print "PAM  search %s (%d)" % (pam,cur_dist)
#    
    return(best_guide_SNP_position, best_pam, all_pam_SNP_positions, all_pam_per_position)


###############################################################################################################
def extract_guide_for_position(snp_pos, chrom, guide_length, edit_max_distance_from_PAM5prime, pam_seqs, downstream_or_upstream, fasta_sequences):
  # snp_position is zero based
  # chrom - should match the fasta sequences chromosomes format
  # guide_length - 20bp (int)
  # edit_max_distance_from_PAM5prime - positive int. 18
  # pam_seqs - list of PAM sequences ['GA]
  # downstream_or_upstream ['downstream'|'upstream']


  # alt_allele - alternative allele  not used in the current versions (waas used to avoid GG -> GA edits)
  alt_allele = "N"
  
  if (downstream_or_upstream == 'downstream'):
    
    SNP_sense = fasta_sequences[chrom][(snp_pos-1):(snp_pos+edit_max_distance_from_PAM5prime+3)]

    pos_in_guide, min_dist_pam,  all_pam_SNP_positions, all_pam_per_position = \
                                  get_best_guide_SNP_position(SNP_sense.seq, pam_seqs, alt_allele)


    if (np.isnan(pos_in_guide)):
      guide = Seq("",generic_dna)
      guide_revcomp = Seq("",generic_dna)
      guide_noPAM = Seq("",generic_dna)
      guide_PAM_p7 = Seq("",generic_dna)
      guide_PAM_m4p3 = Seq("",generic_dna)
    else:
      guide_end_pos = int(snp_pos+2-pos_in_guide)

      # the guide contains the PAM
      guide = fasta_sequences[chrom].seq[(guide_end_pos-guide_length+1-3):(guide_end_pos+1)]
      guide_noPAM = fasta_sequences[chrom].seq[(guide_end_pos-guide_length+1-3):(guide_end_pos+1)-3]
      guide_PAM_p7 = fasta_sequences[chrom].seq[(guide_end_pos-guide_length+1-3):(guide_end_pos+1)+7]
      guide_PAM_m4p3 = fasta_sequences[chrom].seq[(guide_end_pos-guide_length+1-3-4):(guide_end_pos+1)+3]

      guide_revcomp = guide.reverse_complement()

  elif (downstream_or_upstream == 'upstream'):
    
    SNP_antisense = fasta_sequences[chrom][(snp_pos-edit_max_distance_from_PAM5prime-3):(snp_pos+2)]
    SNP_antisense = SNP_antisense.reverse_complement()

    # TODO - is OK for alternative longer than one
    #alt_revcomp = Seq(str(record.ALT[0]),generic_dna).reverse_complement()
    alt_revcomp = alt_allele

    pos_in_guide, min_dist_pam,  all_pam_SNP_positions, all_pam_per_position = \
    get_best_guide_SNP_position(SNP_antisense.seq, pam_seqs, alt_revcomp)

    if (np.isnan(pos_in_guide)):
      guide = Seq("",generic_dna)
      guide_revcomp = Seq("",generic_dna)
      guide_noPAM = Seq("",generic_dna)
      guide_PAM_p7 = Seq("",generic_dna)
      guide_PAM_m4p3 = Seq("",generic_dna)
    else:
      guide_strart_pos = int(snp_pos-2+pos_in_guide)

      # the guide contains the PAM
      guide = fasta_sequences[chrom].seq[(guide_strart_pos):(guide_strart_pos+guide_length+3)]
      #print guide
      guide_revcomp = guide.reverse_complement()
      #print guide_revcomp

      guide_noPAM = (fasta_sequences[chrom].seq[(guide_strart_pos+3):(guide_strart_pos+guide_length+3)]).reverse_complement()
      guide_PAM_p7 =  (fasta_sequences[chrom].seq[(guide_strart_pos-7):(guide_strart_pos+guide_length+3)]).reverse_complement()
      guide_PAM_m4p3 =  (fasta_sequences[chrom].seq[(guide_strart_pos-3):(guide_strart_pos+guide_length+3+4)]).reverse_complement()

  else:
    print(extract_guide_for_position  + "downstream_or_upstream can only be downstream or upstream")
    raise ValueError(extract_guide_for_position  + "downstream_or_upstream can only be downstream or upstream")
  
  
  if np.isnan(pos_in_guide):
    pos_in_guide_2print = None
  else:
    pos_in_guide_2print = int(pos_in_guide)
    
    
  return(  {'guide' :  guide, \
            'guide_revcomp' : guide_revcomp, \
            'guide_noPAM' : guide_noPAM, \
            'guide_PAM_p7' : guide_PAM_p7, \
            'guide_PAM_m4p3' : guide_PAM_m4p3, \
            'pos_in_guide_2print' : pos_in_guide_2print, \
            'pos_in_guide' : pos_in_guide, \
            'min_dist_pam' : min_dist_pam, \
            'all_pam_SNP_positions' : all_pam_SNP_positions, \
            'all_pam_per_position' : all_pam_per_position} )


###############################################################################################################
def extract_guides_for_positions(input_pos_tab_filename, genome_fasta_file, \
                                 output_filename_SNP, output_filename_GUIDE, \
                                 guide_length, pam_seqs, edit_max_distance_from_PAM5prime):
    # loading genome fasta file
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(genome_fasta_file),'fasta'))
    # set genome to uppercase
    fasta_sequences = {chrom : sequence.upper() for chrom, sequence in fasta_sequences.items()}
    
    input_position_forGuides_df = pd.read_table(input_pos_tab_filename, sep='\t', na_values = "None")

    
    with open(output_filename_SNP, mode='w') as out_tab_file:
        
        out_tab_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                           ("SNPID","CHROM","POS","REF","ALT", "ALT_revcomp","Guide_exists","num_possible_guides",
                            "pos_in_guide_downstream", "min_dist_pam_downstream", "guide_downstream", "guide_downstream_revcomp",
                            "guide_downstream_noPAM", "guide_downstream_PAM_p7","guide_downstream_PAM_m4p3",
                            "downstream_all_pam_SNP_positions", "downstream_all_pam_per_position",
                            "pos_in_guide_upstream", "min_dist_pam_upstream", "guide_upstream", "guide_upstream_revcomp",
                            "guide_upstream_noPAM", "guide_upstream_PAM_p7", "guide_upstream_PAM_m4p3",
                            "upstream_all_pam_SNP_positions", "upstream_all_pam_per_position",
                            "PlusMinus25") ) 


        with open(output_filename_GUIDE, mode='w') as out_tab_file_guides:

            out_tab_file_guides.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                               ("SNPID","CHROM","POS","REF","ALT", 
                                "GUIDE_ID", "index_in_snp", "upstream_or_downstream", "SNP_pos_in_guide",
                                "guide_noPAM", "guide_PAM_7bp", "guide_PAM_m4p3bp") )

            guide_id = 0
            found_pam=0
            for i, row in input_position_forGuides_df.iterrows():
                
                #print row['POS'], row['CHROM']
                
                snp_pos = int(row['POS']) - 1
                chrom = str(row['CHROM'])
                
                #print str(snp_pos)
                #print chrom
                #print type(snp_pos)
                #print type(chrom)

                if (i % 1000 == 0):
                  sys.stderr.write("Parsing line: %i, found: %d, #guides: %d\n" % (i, found_pam, guide_id))


                extracted_guides_downstream = extract_guide_for_position(snp_pos, chrom, guide_length, \
                                                              edit_max_distance_from_PAM5prime, pam_seqs, 'downstream', fasta_sequences)
                
                extracted_guides_upstream = extract_guide_for_position(snp_pos, chrom, guide_length, \
                                                              edit_max_distance_from_PAM5prime, pam_seqs, 'upstream', fasta_sequences)

                
                ################################################################################
                # writing to the SNP file
                ################################################################################
                

                Guide_exists = False
                if (not np.isnan(extracted_guides_downstream['pos_in_guide']) or not np.isnan(extracted_guides_upstream['pos_in_guide'])):
                    found_pam=found_pam+1
                    Guide_exists = True


                num_possible_guides = len(extracted_guides_downstream['all_pam_SNP_positions']) + \
                                      len(extracted_guides_upstream['all_pam_SNP_positions'])

                out_tab_file.write("%d\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                                   (i,str(row['CHROM']),int(row['POS']),
                                    fasta_sequences[chrom].seq[snp_pos:(snp_pos+1)], fasta_sequences[chrom].seq[snp_pos:(snp_pos+1)], fasta_sequences[chrom].seq[snp_pos:(snp_pos+1)],
                                    Guide_exists, num_possible_guides,
                                    
                                    str(extracted_guides_downstream['pos_in_guide_2print']), extracted_guides_downstream['min_dist_pam'], 
                                    extracted_guides_downstream['guide'], extracted_guides_downstream['guide_revcomp'],
                                    extracted_guides_downstream['guide_noPAM'], extracted_guides_downstream['guide_PAM_p7'],extracted_guides_downstream['guide_PAM_m4p3'],
                                    str(extracted_guides_downstream['all_pam_SNP_positions'].astype(int)), str(extracted_guides_downstream['all_pam_per_position']),
                                    
                                    str(extracted_guides_upstream['pos_in_guide_2print']), extracted_guides_upstream['min_dist_pam'], 
                                    extracted_guides_upstream['guide'], extracted_guides_upstream['guide_revcomp'],
                                    extracted_guides_upstream['guide_noPAM'], extracted_guides_upstream['guide_PAM_p7'],extracted_guides_upstream['guide_PAM_m4p3'],
                                    str(extracted_guides_upstream['all_pam_SNP_positions'].astype(int)), str(extracted_guides_upstream['all_pam_per_position']),                                    
                                    fasta_sequences[chrom].seq[(snp_pos-25):(snp_pos+25)]) )
                           
                ################################################################################
                # writing to guides files
                ################################################################################
 
                downstream_all_pam_SNP_positions = extracted_guides_downstream['all_pam_SNP_positions']
                upstream_all_pam_SNP_positions = extracted_guides_upstream['all_pam_SNP_positions']
                

                index_in_snp = 0

                for cur_SNP_pos_in_guide in downstream_all_pam_SNP_positions.astype(int):
                    guide_id = guide_id + 1
                    index_in_snp = index_in_snp + 1
                    cur_upstream_or_downstream = "downstream"
                    cur_guide_downstream_end_pos = int(snp_pos+2-cur_SNP_pos_in_guide)

                    cur_guide_noPAM =  fasta_sequences[chrom].seq[(cur_guide_downstream_end_pos-guide_length+1-3):(cur_guide_downstream_end_pos+1)-3]
                    cur_guide_PAM_p7 = fasta_sequences[chrom].seq[(cur_guide_downstream_end_pos-guide_length+1-3):(cur_guide_downstream_end_pos+1)+7]
                    cur_guide_PAM_m4p3 = fasta_sequences[chrom].seq[(cur_guide_downstream_end_pos-guide_length+1-3-4):(cur_guide_downstream_end_pos+1)+3]

                    out_tab_file_guides.write("%d\t%s\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n" % 
                                              (i,str(row['CHROM']),int(row['POS']),
                                               fasta_sequences[chrom].seq[int(row['POS']):(int(row['POS'])+1)],fasta_sequences[chrom].seq[int(row['POS']):(int(row['POS'])+1)],
                                               guide_id, index_in_snp, cur_upstream_or_downstream, cur_SNP_pos_in_guide,
                                               cur_guide_noPAM, cur_guide_PAM_p7,cur_guide_PAM_m4p3) )
                    #guide_list.append(guide_id)
                    #guide_id_list_p7.append(str(cur_guide_PAM_p7))
                    #guide_id_list_m4p3.append(str(cur_guide_PAM_m4p3))


                for cur_SNP_pos_in_guide in upstream_all_pam_SNP_positions.astype(int):
                    guide_id = guide_id + 1
                    index_in_snp = index_in_snp + 1
                    cur_upstream_or_downstream = "upstream"
                    cur_guide_upstream_strart_pos = int(snp_pos-2+cur_SNP_pos_in_guide)

                    cur_guide_noPAM =   (fasta_sequences[chrom].seq[(cur_guide_upstream_strart_pos+3):(cur_guide_upstream_strart_pos+guide_length+3)]).reverse_complement()
                    cur_guide_PAM_p7 =  (fasta_sequences[chrom].seq[(cur_guide_upstream_strart_pos-7):(cur_guide_upstream_strart_pos+guide_length+3)]).reverse_complement()
                    cur_guide_PAM_m4p3 =  (fasta_sequences[chrom].seq[(cur_guide_upstream_strart_pos-3):(cur_guide_upstream_strart_pos+guide_length+3+4)]).reverse_complement()

                    out_tab_file_guides.write("%d\t%s\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n" % 
                                              (i,str(row['CHROM']),int(row['POS']),
                                               fasta_sequences[chrom].seq[int(row['POS']):(int(row['POS'])+1)],fasta_sequences[chrom].seq[int(row['POS']):(int(row['POS'])+1)],
                                               guide_id, index_in_snp, cur_upstream_or_downstream, cur_SNP_pos_in_guide,
                                               cur_guide_noPAM, cur_guide_PAM_p7, cur_guide_PAM_m4p3) )
                    #guide_list.append(guide_id)
                    #guide_id_list_p7.append(str(cur_guide_PAM_p7))
                    #guide_id_list_m4p3.append(str(cur_guide_PAM_m4p3))



                #if (i>1):
                #    break

            sys.stderr.write("Finish parsing input file: %i, found: %d, #guides: %d\n" % (i, found_pam, guide_id))


###############################################################################################################
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
    guides_df['guide_map_mismatch_0'] = 0
    guides_df['guide_map_mismatch_1'] = 0
    guides_df['guide_map_mismatch_2'] = 0
    guides_df['guide_map_mismatch_3'] = 0
    
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
                raise ValueError("off target mapping: bowtie2 line does not contain XM or NM:" + line)
            
            
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
def cal_guide_features(input_guides_table_filename, input_genome_fasta_filename, output_guides_with_features_table_filename = "", 
                       PAM_seq = 'GA', excluded_seqs = ['TTTTT'], donor_length = 100 ,BOWTIE_exe = which('bowtie2'),
                        run_azimuth = False):
    """
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

    """

    
    # fixing the mapping command to be the default bowtie
    if not BOWTIE_exe or len(BOWTIE_exe) < 3:
        
        print("input bowtie2 is null, using which function to find path")
        BOWTIE_exe = which('bowtie2')
        
        if not BOWTIE_exe:
            raise ValueError("bowtie2 is unknown (which returns None, make sure unix >>which bowtie2)")
        print('###')
        


    # defining the bowtie mapping cmd
    bowtie_cmd = shlex.split(BOWTIE_exe + ' -x ' + os.path.splitext(input_genome_fasta_filename)[0] + ' -U - -f -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 --gbar 3 --end-to-end -k 30 --no-head -t  --rdg 10,6 --rfg 10,6')
    mapping_cmd = bowtie_cmd
    
    # TODO remove this printing 
    print("Using this command for mapping off targets: " + BOWTIE_exe + ' -x ' + os.path.splitext(input_genome_fasta_filename)[0] + ' -U - -f -D 20 -R 3 -N 1 -L 10 -i S,1,0.50 --gbar 3 --end-to-end -k 30 --no-head -t  --rdg 10,6 --rfg 10,6')

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
    out_guides_with_features_df = eval_guides_off_targets(out_guides_with_features_df, 
                                                          genome_seq, guide_id_col_name = 'guide_id', guide_col_name = 'guide_noPAM', 
                                                          pam_2nt_seq = PAM_seq, mapping_cmd = mapping_cmd,
                                                          input_genome_fasta_filename = input_genome_fasta_filename)
    
    
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


###############################################################################################################
def extract_guides_for_snps(input_snps_vcf_filename, input_genome_fasta_filename, 
                            output_SNP_table_filename, output_guides_table_filename,
                            pam_seqs = ["GG"], guide_length = 20, edit_max_distance_from_PAM5prime = 18,
                            var_id_prefix = ""):
    """
    extract guides for SNPs.
    Write SNP table and guides table 
    
    """
    # loading the VCF SNPs file
    vcf_reader = vcf.Reader(open(input_snps_vcf_filename, 'r'))
    
    # loading genome fasta file
    genome_seq = SeqIO.to_dict(SeqIO.parse(open(input_genome_fasta_filename),'fasta', alphabet=generic_dna))
    # set genome to uppercase
    genome_seq = {chrom : sequence.upper() for chrom, sequence in genome_seq.items()}
 

    print('---------------------- extracting guides for SNPs -------------------------------')

    with open(output_SNP_table_filename, mode='w') as out_tab_file_SNPs:

        # SNP table header
        out_tab_file_SNPs.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                           ("var_id","chrom","SNP_chr_pos","REF","ALT", "ALT_revcomp","Guide_exists","num_possible_guides",
                            "pos_in_guide_downstream", "min_dist_pam_downstream", "guide_downstream", "guide_downstream_revcomp",
                            "guide_downstream_noPAM", "guide_downstream_PAM_p7","guide_downstream_PAM_m4p3",
                            "downstream_all_pam_SNP_positions", "downstream_all_pam_per_position",
                            "pos_in_guide_upstream", "min_dist_pam_upstream", "guide_upstream", "guide_upstream_revcomp",
                            "guide_upstream_noPAM", "guide_upstream_PAM_p7", "guide_upstream_PAM_m4p3",
                            "upstream_all_pam_SNP_positions", "upstream_all_pam_per_position",
                            "PlusMinus25") ) 


        with open(output_guides_table_filename, mode='w') as out_tab_file_guides:

            # guide table header
            out_tab_file_guides.write( 
                                ("%s\t%s\t%s\t%s\t%s\t" +  
                                 "%s\t%s\t%s\t%s\t%s\t" +
                                 "%s\t%s\t" +
                                 "%s\t%s\t%s\t%s\t\n") % 
                               ("var_id","chrom","SNP_chr_pos","REF","ALT", 
                                "guide_id", "guide_strand", "guide0_chr_pos", "guide_cut_chr_pos", "SNP_pos_in_guide",
                                "guide_index_in_snp", "upstream_or_downstream", 
                                "guide", "guide_noPAM", "guide_PAM_7bp", "guide_PAM_m4p3bp") )

            # iterating over the SNPs
            i=0
            guide_id = 0
            found_pam=0
            for record in vcf_reader:
                i=i+1

                if (i % 5000 == 0):
                  print("Parsing SNP file line: %i, found guides for: %d, #total guides: %d" % (i, found_pam, guide_id))

                snp_pos = record.POS-1

                ################################################################################

                #####################
                # downstream guide
                #####################
                SNP_sense_downstream = genome_seq[record.CHROM][(snp_pos-1):(snp_pos+edit_max_distance_from_PAM5prime+3)]

                pos_in_guide_downstream, min_dist_pam_downstream,  downstream_all_pam_SNP_positions, downstream_all_pam_per_position = \
                get_best_guide_SNP_position(SNP_sense_downstream.seq, pam_seqs, record.ALT)


                if (np.isnan(pos_in_guide_downstream)):
                    guide_downstream = Seq("",generic_dna)
                    guide_downstream_revcomp = Seq("",generic_dna)
                    guide_downstream_noPAM = Seq("",generic_dna)
                    guide_downstream_PAM_p7 = Seq("",generic_dna)
                    guide_downstream_PAM_m4p3 = Seq("",generic_dna)
                else:
                    guide_downstream_end_pos = int(snp_pos+2-pos_in_guide_downstream)

                    # the guide contains the PAM
                    guide_downstream          = genome_seq[record.CHROM].seq[(guide_downstream_end_pos-guide_length+1-3):(guide_downstream_end_pos+1)]
                    guide_downstream_noPAM    = genome_seq[record.CHROM].seq[(guide_downstream_end_pos-guide_length+1-3):(guide_downstream_end_pos+1)-3]
                    guide_downstream_PAM_p7   = genome_seq[record.CHROM].seq[(guide_downstream_end_pos-guide_length+1-3):(guide_downstream_end_pos+1)+7]
                    guide_downstream_PAM_m4p3 = genome_seq[record.CHROM].seq[(guide_downstream_end_pos-guide_length+1-3-4):(guide_downstream_end_pos+1)+3]

                    guide_downstream_revcomp = guide_downstream.reverse_complement()


                #####################
                # upstream guide
                #####################

                #print "******"

                SNP_antisense_upstream = genome_seq[record.CHROM][(snp_pos-edit_max_distance_from_PAM5prime-2):(snp_pos+2)]
                #print SNP_antisense_upstream.seq
                SNP_antisense_upstream = SNP_antisense_upstream.reverse_complement()
                #print SNP_antisense_upstream.seq

                # TODO - is OK for alternative longer than one
                alt_revcomp = Seq(str(record.ALT[0]),generic_dna).reverse_complement()

                pos_in_guide_upstream, min_dist_pam_upstream,  upstream_all_pam_SNP_positions, upstream_all_pam_per_position = \
                get_best_guide_SNP_position(SNP_antisense_upstream.seq, pam_seqs, alt_revcomp)

                if (np.isnan(pos_in_guide_upstream)):
                    guide_upstream = Seq("",generic_dna)
                    guide_upstream_revcomp = Seq("",generic_dna)
                    guide_upstream_noPAM = Seq("",generic_dna)
                    guide_upstream_PAM_p7 = Seq("",generic_dna)
                    guide_upstream_PAM_m4p3 = Seq("",generic_dna)
                else:
                    guide_upstream_strart_pos = int(snp_pos-2+pos_in_guide_upstream)


                    # the guide contains the PAM
                    guide_upstream = genome_seq[record.CHROM].seq[(guide_upstream_strart_pos):(guide_upstream_strart_pos+guide_length+3)]
                    guide_upstream_revcomp = guide_upstream.reverse_complement()


                    guide_upstream_noPAM = (genome_seq[record.CHROM].seq[(guide_upstream_strart_pos+3):(guide_upstream_strart_pos+guide_length+3)]).reverse_complement()
                    guide_upstream_PAM_p7 =  (genome_seq[record.CHROM].seq[(guide_upstream_strart_pos-7):(guide_upstream_strart_pos+guide_length+3)]).reverse_complement()
                    guide_upstream_PAM_m4p3 =  (genome_seq[record.CHROM].seq[(guide_upstream_strart_pos-3):(guide_upstream_strart_pos+guide_length+3+4)]).reverse_complement()


                ################################################################################

                # printing 

                Guide_exists = False
                if (not np.isnan(pos_in_guide_downstream) or not np.isnan(pos_in_guide_upstream)):
                    found_pam=found_pam+1
                    Guide_exists = True


                if np.isnan(pos_in_guide_downstream):
                    pos_in_guide_downstream_2print = None
                else:
                    pos_in_guide_downstream_2print = int(pos_in_guide_downstream)



                if np.isnan(pos_in_guide_upstream):
                    pos_in_guide_upstream_2print = None
                else:
                    pos_in_guide_upstream_2print = int(pos_in_guide_upstream)

                num_possible_guides = len(downstream_all_pam_SNP_positions) + len(upstream_all_pam_SNP_positions)

                out_tab_file_SNPs.write("%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                                   (write_var_id(record=record,var_id_prefix=var_id_prefix,i=i),record.CHROM,record.POS,record.REF,record.ALT, alt_revcomp,Guide_exists, num_possible_guides,
                                    str(pos_in_guide_downstream_2print), min_dist_pam_downstream, guide_downstream, guide_downstream_revcomp,
                                    guide_downstream_noPAM, guide_downstream_PAM_p7,guide_downstream_PAM_m4p3,
                                    str(downstream_all_pam_SNP_positions.astype(int)), str(downstream_all_pam_per_position),
                                    str(pos_in_guide_upstream_2print), min_dist_pam_upstream, guide_upstream, guide_upstream_revcomp,
                                    guide_upstream_noPAM, guide_upstream_PAM_p7,guide_upstream_PAM_m4p3,
                                    str(upstream_all_pam_SNP_positions.astype(int)), str(upstream_all_pam_per_position),
                                    genome_seq[record.CHROM].seq[(snp_pos-25):(snp_pos+25)]) )


                ######################
                # writing to guides files
                ######################

                index_in_snp = 0

                for cur_SNP_pos_in_guide in downstream_all_pam_SNP_positions.astype(int):
                    guide_id = guide_id + 1
                    index_in_snp = index_in_snp + 1
                    cur_upstream_or_downstream = "downstream"
                    cur_guide_strand = "+"
                    cur_guide_downstream_end_pos = int(snp_pos+2-cur_SNP_pos_in_guide)

                    cur_guide          =  genome_seq[record.CHROM].seq[(cur_guide_downstream_end_pos-guide_length+1-3):(cur_guide_downstream_end_pos+1)]
                    cur_guide_noPAM    =  genome_seq[record.CHROM].seq[(cur_guide_downstream_end_pos-guide_length+1-3):(cur_guide_downstream_end_pos+1)-3]
                    cur_guide_PAM_p7   = genome_seq[record.CHROM].seq[(cur_guide_downstream_end_pos-guide_length+1-3):(cur_guide_downstream_end_pos+1)+7]
                    cur_guide_PAM_m4p3 = genome_seq[record.CHROM].seq[(cur_guide_downstream_end_pos-guide_length+1-3-4):(cur_guide_downstream_end_pos+1)+3]

                    cur_guide0_chr_pos = int(cur_guide_downstream_end_pos-2)
                    # the bp downstream to the cut
                    cur_guide_cut_chr_pos = int((cur_guide_downstream_end_pos-2)-3)
                    
                     

                    out_tab_file_guides.write( ("%s\t%s\t%d\t%s\t%s\t" + 
                                              "%s\t%s\t%d\t%d\t%d\t" +
                                              "%d\t%s\t" +
                                              "%s\t%s\t%s\t%s\n") % 
                                              (write_var_id(record=record,var_id_prefix=var_id_prefix,i=i), record.CHROM, record.POS, record.REF, record.ALT,
                                               ("guide_" + str(guide_id)), cur_guide_strand, cur_guide0_chr_pos, cur_guide_cut_chr_pos, cur_SNP_pos_in_guide,
                                               index_in_snp, cur_upstream_or_downstream,
                                               str(cur_guide), str(cur_guide_noPAM), str(cur_guide_PAM_p7),str(cur_guide_PAM_m4p3) ) )


                for cur_SNP_pos_in_guide in upstream_all_pam_SNP_positions.astype(int):
                    guide_id = guide_id + 1
                    index_in_snp = index_in_snp + 1
                    cur_upstream_or_downstream = "upstream"
                    cur_guide_strand = "-"
                    cur_guide_upstream_strart_pos = int(snp_pos-2+cur_SNP_pos_in_guide)

                    cur_guide          =   (genome_seq[record.CHROM].seq[(cur_guide_upstream_strart_pos):(cur_guide_upstream_strart_pos+guide_length+3)]).reverse_complement()
                    cur_guide_noPAM    =   (genome_seq[record.CHROM].seq[(cur_guide_upstream_strart_pos+3):(cur_guide_upstream_strart_pos+guide_length+3)]).reverse_complement()
                    cur_guide_PAM_p7   =  (genome_seq[record.CHROM].seq[(cur_guide_upstream_strart_pos-7):(cur_guide_upstream_strart_pos+guide_length+3)]).reverse_complement()
                    cur_guide_PAM_m4p3 =  (genome_seq[record.CHROM].seq[(cur_guide_upstream_strart_pos-3):(cur_guide_upstream_strart_pos+guide_length+3+4)]).reverse_complement()

                    cur_guide0_chr_pos = int(cur_guide_upstream_strart_pos+2)
                    # the bp downstream to the cut
                    cur_guide_cut_chr_pos = int((cur_guide_upstream_strart_pos+2)+4)
                    
                    
                    out_tab_file_guides.write( ("%s\t%s\t%d\t%s\t%s\t" + 
                                              "%s\t%s\t%d\t%d\t%d\t" +
                                              "%d\t%s\t" +
                                              "%s\t%s\t%s\t%s\n") % 
                                              (write_var_id(record=record,var_id_prefix=var_id_prefix,i=i), record.CHROM, record.POS, record.REF, record.ALT,
                                               ("guide_" + str(guide_id)), cur_guide_strand, cur_guide0_chr_pos, cur_guide_cut_chr_pos, cur_SNP_pos_in_guide,
                                               index_in_snp, cur_upstream_or_downstream,
                                               str(cur_guide), str(cur_guide_noPAM), str(cur_guide_PAM_p7),str(cur_guide_PAM_m4p3) ) )
                    


                #if (i>1):
                #    break

            print("Finish parsing VCF: %i, found: %d, #guides: %d\n" % (i, found_pam, guide_id))
    
    print("---------------------------- Done extracting guides for SNPs --------------------------")



###############################################################################################################
# functions for donor sequences design for SNP guides
###############################################################################################################

###############################################################################################################
def read_snp_donor_design_table(input_design_table_filename):
    
    design_df = pd.read_table(input_design_table_filename, sep='\t')
    design_df['donor_seq_offsets'] = design_df['donor_seq_offsets'].apply(ast.literal_eval)
    design_df['filter_in'] = design_df['filter_in'].str.strip()
    design_df['filter_out'] = design_df['filter_out'].str.strip()
    design_df['donor_mut_type'] = design_df['donor_mut_type'].str.strip()
    
    return(design_df)

###############################################################################################################
def get_design_filter_guide_ids(design_row, snps_df, guides_df):
    
    
    filter_in_strs = design_row['filter_in']
    filter_out_strs = design_row['filter_out']
    
    filter_guide_ids = guides_df['guide_id'].copy()
    
    
    if  filter_in_strs and type(filter_in_strs) is str and filter_in_strs != "None":
        #print "in: " + filter_in_strs
        for filter_str in filter_in_strs.split(','):
            cur_filter_str = filter_str.strip()
            cur_filter_guide_ids = get_design_filter_str_guide_ids(cur_filter_str, snps_df, guides_df)
            filter_guide_ids = filter_guide_ids[filter_guide_ids.isin(cur_filter_guide_ids)]
    
    if  filter_out_strs and type(filter_out_strs) is str and filter_out_strs != "None":
        #print "out: " + filter_out_strs
        for filter_str in filter_out_strs.split(','):
            cur_filter_str = filter_str.strip()
            cur_filter_guide_ids = get_design_filter_str_guide_ids(cur_filter_str, snps_df, guides_df)
            filter_guide_ids = filter_guide_ids[~filter_guide_ids.isin(cur_filter_guide_ids)]      
    
    return(filter_guide_ids)
    
###############################################################################################################
def get_design_filter_str_guide_ids(filter_str, snps_df, guides_df):
    
    #print "in get_design_filter_str_guide_ids"
    filter_fields = filter_str.split(':')
    
    if filter_fields[0] == "VAR":
        var_prefix_list = filter_fields[1].split(';')
        selection = guides_df['var_id'].apply(lambda x: x.split('_')[0]).isin(var_prefix_list)
        ret_guide_ids = guides_df.loc[selection, 'guide_id']
    
    #TODO - debug the type conversion
    elif filter_fields[0] == "SNP":
        ret_guide_ids = guides_df['guide_id'][ guides_df['var_id'].isin( snps_df['var_id'][snps_df[filter_fields[1]] ==  
                                                convert_str_to_pandas_dtype_comparable(snps_df[filter_fields[1]],filter_fields[2])] ) ]
    elif filter_fields[0] == "GUIDE":
        ret_guide_ids = guides_df['guide_id'][ guides_df[ filter_fields[1]] ==  
                                                convert_str_to_pandas_dtype_comparable(guides_df[filter_fields[1]],filter_fields[2])]
    else:
        raise ValueError("unknown matrix type in donor design filter str:" + filter_str)
    
    return(ret_guide_ids)

###############################################################################################################
def convert_str_to_pandas_dtype_comparable(pd_series, val_str):
    """
    convert a string to a type that will be comparable to the 
    """
    return (get_pandas_col_python_type(pd_series)(val_str))

###############################################################################################################
def get_pandas_col_python_type(pd_series):
    # TODO - add all possible types
    if pd_series.dtype.name in ['int8','int16','int32','int64']:
        return type(int(1))
    elif pd_series.dtype.name in ['float16','float32','float64','float128']:
        return type(float(1.0))
    elif pd_series.dtype.name in ['string','object']:
        return type(str(1.0))
    elif pd_series.dtype.name in ['bool']:
        return type(bool(True))
    else:
        raise ValueError('Unknown dtype: %s' % (pd_series.dtype.name))
###############################################################################################################
def get_donor_mut_for_SNP_guides_old2020(guide_ids, snps_df, guides_df, genome_seq,
                                donor_mut_type, donor_length, excluded_seqs,
                                donor_seq_offsets, min_dist_cut_to_donor_edge,
                                set_name, donor_id):
    
    out_guide_donor_df =  pd.DataFrame(data=None)
    
    for guide_id in guide_ids:
        
        if donor_seq_offsets is None:
            #print "donor offsets is None" #DEBUG
            cur_guide_donor_df, donor_id = get_donor_mut_for_SNP_guide(guide_id, snps_df, guides_df, genome_seq,
                                                               donor_mut_type, donor_length, excluded_seqs, None, min_dist_cut_to_donor_edge,
                                                               set_name, donor_id)
            out_guide_donor_df = out_guide_donor_df.append(cur_guide_donor_df,ignore_index=True)
            
            if (donor_id % 1000 == 0):
                print("Designing donor number: %d" % (donor_id))
        else:
            for donor_seq_offset in donor_seq_offsets:
                cur_guide_donor_df, donor_id = get_donor_mut_for_SNP_guide(guide_id, snps_df, guides_df, genome_seq,
                                                               donor_mut_type, donor_length, excluded_seqs, donor_seq_offset, min_dist_cut_to_donor_edge,
                                                               set_name, donor_id)
                out_guide_donor_df = out_guide_donor_df.append(cur_guide_donor_df,ignore_index=True)
                if (donor_id % 1000 == 0):
                    print("Designing donor number: %d" % (donor_id))
    
    return (out_guide_donor_df, donor_id)

###############################################################################################################
def get_donor_mut_for_SNP_guides(guide_ids, snps_df, guides_df, genome_seq,
                                donor_mut_type, donor_length, excluded_seqs,
                                donor_seq_offsets, min_dist_5prime_arm, min_dist_3prime_arm,
                                set_name, donor_id):
    """
    modified to accommodate different asymmetrical minimum distances from cut to donor edge
    """
    
    out_guide_donor_df =  pd.DataFrame(data=None)
    
    for guide_id in guide_ids:
        
        if donor_seq_offsets is None:
            #print "donor offsets is None" #DEBUG
            cur_guide_donor_df, donor_id = get_donor_mut_for_SNP_guide(guide_id, snps_df, guides_df, genome_seq,
                                                               donor_mut_type, donor_length, excluded_seqs, None, min_dist_5prime_arm, min_dist_3prime_arm, set_name, donor_id)
            out_guide_donor_df = out_guide_donor_df.append(cur_guide_donor_df,ignore_index=True)
            
            if (donor_id % 1000 == 0):
                print("Designing donor number: %d" % (donor_id))
        else:
            for donor_seq_offset in donor_seq_offsets:
                cur_guide_donor_df, donor_id = get_donor_mut_for_SNP_guide(guide_id, snps_df, guides_df, genome_seq,
                                                               donor_mut_type, donor_length, excluded_seqs, donor_seq_offset, min_dist_5prime_arm, min_dist_3prime_arm, set_name, donor_id)
                out_guide_donor_df = out_guide_donor_df.append(cur_guide_donor_df,ignore_index=True)
                if (donor_id % 1000 == 0):
                    print("Designing donor number: %d" % (donor_id))
    
    return (out_guide_donor_df, donor_id)

###############################################################################################################
def get_donor_mut_for_SNP_guide_old2_1_2018(guide_id, snps_df, guides_df, genome_seq,
                                donor_mut_type, donor_length, excluded_seqs,
                                donor_seq_offset, min_dist_cut_to_donor_edge,
                                set_name, donor_id):
    """
    donor_mut_type - REF2ALT
    """
    out_guide_donor_df = pd.DataFrame(data=None)
    
    
    # single row of the guide
    guides_df = guides_df[guides_df['guide_id'] == guide_id]

    # general
    guide_cut_chr_pos = int(guides_df['guide_cut_chr_pos'].iloc[0])
    guide0_chr_pos = int(guides_df['guide0_chr_pos'].iloc[0])
    guide_is_negative_strand = (guides_df['guide_strand'].iloc[0] == '-')
    guide_chrom = str(guides_df['chrom'].iloc[0])
    
    # specific to SNPs
    var_id = str(guides_df['var_id'].iloc[0])
    SNP_pos_in_guide = int(guides_df['SNP_pos_in_guide'].iloc[0])
    snp_ref = str(guides_df['REF'].iloc[0])
    snp_alt = str((guides_df['ALT'].iloc[0])) # assuming '[]' was removed from the string
    snp_chr_pos = int(guides_df['SNP_chr_pos'].iloc[0]) - 1 # TODO this is one based? (VCF)
    
    
    donor_mut_type_splt = donor_mut_type.split('_')
    donor_mut_name = donor_mut_type_splt[0]
    


    donor_seq = Seq("", generic_dna)
    
    
    if  (donor_mut_name == 'REF2ALT'):
        
        #len_diff_alt_ref = len(snp_alt) - len(snp_ref)
        shift_dis_from_cut =  int(snp_chr_pos - guide_cut_chr_pos + 1) # + 1 * (!guide_is_negative_strand)  )
  
        
        # None means try to find a donor sequence with no excluded seqs
        if donor_seq_offset is None:
            
            #DEBUG
            #print "Looking for no excluding for guide_id: " + guide_id
            #print "guide_is_negative_strand"
            #print guide_is_negative_strand
            
            donor_seq_offset  = 0
            donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
            # start from looking at the left part
            left_seq = str(genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ]) 
            
            #print "1 donor_seq_offset:  %d" % (donor_seq_offset)  #DEBUG
            #print "donor_length %d" % (donor_length)
            #print "len(snp_alt) %d" % (len(snp_alt))
            #print "shift_dis_from_cut %d" % (shift_dis_from_cut)
            #print "donor_nt_add_left %d" % (donor_nt_add_left)
            #print left_seq
            
            # (may still leave excluded seq overlapping the SNP) - now it is in the middle
            cur_excluded_seq_end = 0
            for m in re.finditer('|'.join(excluded_seqs),left_seq + str(Seq(snp_alt,generic_dna))  ):
                #print m.end() 
                if ( m.end() < len(left_seq)) and m.end() > cur_excluded_seq_end:
                    cur_excluded_seq_end = int( np.floor( (m.end()+m.start())/2 ) )
                    
            #print "cur_excluded_seq_end:" + str(cur_excluded_seq_end)
            
            # removing the part with the excluded seqs (TODO - can do better, remove only part of the excluded_seq)
            if cur_excluded_seq_end>0:
                #print "moving to the right...."
                donor_seq_offset = donor_seq_offset + cur_excluded_seq_end
                donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
                left_seq = str(genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ])
                # after moving left no problem on the right side
                cur_excluded_seq_end = 0

            #print "2 donor_seq_offset:  %d" % (donor_seq_offset)
            
            donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq))
            right_seq = str( genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)])
            
            #print "(snp_chr_pos + len(snp_ref) + donor_nt_add_right) %d" % ((snp_chr_pos + len(snp_ref) + donor_nt_add_right)) 
            
            
            #print "len(left_seq):" + str(len(left_seq))
            #print "len(snp_alt) %d" % (len(snp_alt))
            #print "donor_nt_add_left:" + str(donor_nt_add_left)
            #print left_seq
            #print "donor_nt_add_right:" + str(donor_nt_add_right)
            #print right_seq
            #print "len(right_seq):" + str(len(right_seq))
            
            
            cur_excluded_seq_strat = 0
            for m in re.finditer('|'.join(excluded_seqs),str(Seq(snp_alt,generic_dna)) + right_seq):
                #print m.start()
                cur_suggested_shift = len(right_seq) - int( np.floor( (m.end()+m.start())/2 ) )
                if ( m.start() >= len(snp_alt)  and cur_suggested_shift > cur_excluded_seq_strat ):
                    cur_excluded_seq_strat =  cur_suggested_shift
            
            #print "cur_excluded_seq_strat:" + str(cur_excluded_seq_strat)
            
            # in case the end of the chromosome was reached
            if donor_nt_add_right > len(right_seq):
                print("shifting donor to avoid the end of the chromosome")
                cur_excluded_seq_strat = cur_excluded_seq_strat + (donor_nt_add_right-len(right_seq))
            
            
            #print "cur_excluded_seq_strat:" + str(cur_excluded_seq_strat)
            
            
            if cur_excluded_seq_strat>0:
                #print "moving to the left...."
                donor_seq_offset = donor_seq_offset - cur_excluded_seq_strat
                donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
                left_seq = str(genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ])
                donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq))
                right_seq = str( genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)])
        
                
                #print "2 donor_seq_offset:  %d" % (donor_seq_offset)
                
                # after moving left no problem on the right side
                cur_excluded_seq_strat = 0
            

                cur_excluded_seq_end = 0
                for m in re.finditer('|'.join(excluded_seqs),left_seq + str(Seq(snp_alt,generic_dna))  ):
                    #print m.end() 
                    if ( m.end() < len(left_seq)) and m.end() > cur_excluded_seq_end:
                        cur_excluded_seq_end = int( np.floor( (m.end()+m.start())/2 ) )
                       
                #print "cur_excluded_seq_end:" + str(cur_excluded_seq_end)
   
            
            if (cur_excluded_seq_strat == 0 and cur_excluded_seq_end == 0 and
                len(left_seq) >= min_dist_cut_to_donor_edge and len(right_seq) >= min_dist_cut_to_donor_edge):
                # found ok sequence
                donor_nt_add_left = len(left_seq)
                donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq)) # len(right_seq)
            else:
                warnings.warn("Can NOT get donor with no excluded sequences for guide_id: %s, length left: %d, length right %d (%s,%s) (%d,%d)" % (guide_id, len(left_seq), len(right_seq), left_seq, right_seq, cur_excluded_seq_strat, cur_excluded_seq_end) )
                return (out_guide_donor_df, donor_id)

            #print "found donor without excluded"
            #print int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut + 0 )
            #print donor_nt_add_left
            
            #print "3 donor_seq_offset:  %d" % (donor_seq_offset)
            #donor_seq_offset = int( np.floor((donor_length-len_diff_alt_ref)/2) + shift_dis_from_cut + 0 ) - donor_nt_add_left
            
        else:
            # num of nt to add left and right 
            donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
            donor_nt_add_right = int( (donor_length-len(snp_alt)) - donor_nt_add_left )
        
        
        #print "++++++"
        #print "snp_chr_pos %d" % (snp_chr_pos)
        #print "len(snp_ref) %d" % (len(snp_ref))
        #print "donor_nt_add_right %d" % (donor_nt_add_right)
        
        #print "DEBUG coor %d,  %d,  %d" % ( (snp_chr_pos - donor_nt_add_left) , (snp_chr_pos + len(snp_ref) + donor_nt_add_right), len(genome_seq[guide_chrom].seq))
        #print "DEBUG %s , %s ,  %s" % (str (genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ]) ,
        #                               str(Seq(snp_alt,generic_dna)),
        #                               str(genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)]))
        
        
        # if the SNP is at the edge of the chromosome give warning and do not return a donor        
        if ( (snp_chr_pos - donor_nt_add_left)  <  0 or (snp_chr_pos + len(snp_ref) + donor_nt_add_right) > len(genome_seq[guide_chrom].seq) ):
            warnings.warn("Can NOT get donor for SNP because on the edge of the chromosome. guide_id: %s, length left: %d, length right %d (%s,%s), minimal cooridnate: %d, maximal coordinate: %d, chromosome length: %d" % \
                          (guide_id, len(left_seq), len(right_seq), left_seq, right_seq, (snp_chr_pos - donor_nt_add_left), (snp_chr_pos + len(snp_ref) + donor_nt_add_right) , len(genome_seq[guide_chrom].seq) ) )
            return (out_guide_donor_df, donor_id)
            

            
        donor_seq = genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ] + \
                    Seq(snp_alt,generic_dna)  + \
                    genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)]

        #print "donor_seq before revcomp:  %s" % (str(donor_seq))
        # donor orientation should match the guide
        if (guide_is_negative_strand):
                donor_seq = donor_seq.reverse_complement()
                #print "revcomp the donor seq"
        
        #print "donor_seq after revcomp:  %s" % (str(donor_seq))
        
        
        
        # appending to donor sequences matrix
        donor_id += 1 
        
        cur_donor_line = pd.DataFrame({
            'var_id' : pd.Series(var_id), 'guide_id' : pd.Series(guide_id),
            'donor_id' : pd.Series(guide_id + ':' + donor_mut_type + ':offset' + str(donor_seq_offset) + ':donorID' + str(donor_id) + ':EditPosInGuide' +  str(SNP_pos_in_guide)), 
            'donor_seq': pd.Series(str(donor_seq)), 
            'donor_seq_shift' : pd.Series(int(donor_seq_offset)), 
            'donor_mut_pos_in_guide' : pd.Series(int(SNP_pos_in_guide)), 
            'donor_info_str' : pd.Series("ALT2REF:" + str(snp_ref) + ">" + str(snp_alt)),
            'set_name' : pd.Series(str(set_name))   })
        out_guide_donor_df = out_guide_donor_df.append(cur_donor_line)
        
        
        #print "returning donor"
    else:
        raise ValueError('get_donor_mut_for_SNP_guide unknown donor_mut_type:' + donor_mut_type +  " and donor_mut_name:" + donor_mut_name)
          

    return (out_guide_donor_df, donor_id)


###############################################################################################################
def get_donor_mut_for_SNP_guide_old2020(guide_id, snps_df, guides_df, genome_seq,
                                donor_mut_type, donor_length, excluded_seqs,
                                donor_seq_offset, min_dist_cut_to_donor_edge,
                                set_name, donor_id,
                                to_modify_var_ids = None):
    """
    donor_mut_type -
               REF2ALT,
               towREF2ALT - requires second_var_ids column in the SNPs table,
               multiREF2ALT - reuires a list of SNPs to modify - to_modify_var_ids
    
    
    """
    out_guide_donor_df = pd.DataFrame(data=None)
    
    
    # single row of the guide
    guides_df = guides_df[guides_df['guide_id'] == guide_id]

    # general
    guide_cut_chr_pos = int(guides_df['guide_cut_chr_pos'].iloc[0])
    guide0_chr_pos = int(guides_df['guide0_chr_pos'].iloc[0])
    guide_is_negative_strand = (guides_df['guide_strand'].iloc[0] == '-')
    guide_chrom = str(guides_df['chrom'].iloc[0])
    
    # specific to SNPs
    var_id = str(guides_df['var_id'].iloc[0])
    SNP_pos_in_guide = int(guides_df['SNP_pos_in_guide'].iloc[0])
    snp_ref = str(guides_df['REF'].iloc[0])
    snp_alt = str((guides_df['ALT'].iloc[0])) # assuming '[]' was removed from the string
    snp_chr_pos = int(guides_df['SNP_chr_pos'].iloc[0]) - 1 # TODO this is one based? (VCF)
    
    
    donor_mut_type_splt = donor_mut_type.split('_')
    donor_mut_name = donor_mut_type_splt[0]
    


    donor_seq = Seq("", generic_dna)
    
    
    if  (donor_mut_name == 'REF2ALT'):
        
        #print("designing REF2ALT")
        
        #len_diff_alt_ref = len(snp_alt) - len(snp_ref)
        shift_dis_from_cut =  int(snp_chr_pos - guide_cut_chr_pos + 1) # + 1 * (!guide_is_negative_strand)  )
  
        
        # None means try to find a donor sequence with no excluded seqs
        if donor_seq_offset is None:
            
            #DEBUG
            #print "Looking for no excluding for guide_id: " + guide_id
            #print "guide_is_negative_strand"
            #print guide_is_negative_strand
            
            donor_seq_offset  = 0
            donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
            # start from looking at the left part
            left_seq = str(genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ]) 
            
            #print "1 donor_seq_offset:  %d" % (donor_seq_offset)  #DEBUG
            #print "donor_length %d" % (donor_length)
            #print "len(snp_alt) %d" % (len(snp_alt))
            #print "shift_dis_from_cut %d" % (shift_dis_from_cut)
            #print "donor_nt_add_left %d" % (donor_nt_add_left)
            #print left_seq
            
            # (may still leave excluded seq overlapping the SNP) - now it is in the middle
            cur_excluded_seq_end = 0
            for m in re.finditer('|'.join(excluded_seqs),left_seq + str(Seq(snp_alt,generic_dna))  ):
                #print m.end() 
                if ( m.end() < len(left_seq)) and m.end() > cur_excluded_seq_end:
                    cur_excluded_seq_end = int( np.floor( (m.end()+m.start())/2 ) )
                    
            #print "cur_excluded_seq_end:" + str(cur_excluded_seq_end)
            
            # removing the part with the excluded seqs (TODO - can do better, remove only part of the excluded_seq)
            if cur_excluded_seq_end>0:
                #print "moving to the right...."
                donor_seq_offset = donor_seq_offset + cur_excluded_seq_end
                donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
                left_seq = str(genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ])
                # after moving left no problem on the right side
                cur_excluded_seq_end = 0

            #print "2 donor_seq_offset:  %d" % (donor_seq_offset)
            
            donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq))
            right_seq = str( genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)])
            
            #print "(snp_chr_pos + len(snp_ref) + donor_nt_add_right) %d" % ((snp_chr_pos + len(snp_ref) + donor_nt_add_right)) 
            
            
            #print "len(left_seq):" + str(len(left_seq))
            #print "len(snp_alt) %d" % (len(snp_alt))
            #print "donor_nt_add_left:" + str(donor_nt_add_left)
            #print left_seq
            #print "donor_nt_add_right:" + str(donor_nt_add_right)
            #print right_seq
            #print "len(right_seq):" + str(len(right_seq))
            
            
            cur_excluded_seq_strat = 0
            for m in re.finditer('|'.join(excluded_seqs),str(Seq(snp_alt,generic_dna)) + right_seq):
                #print m.start()
                cur_suggested_shift = len(right_seq) - int( np.floor( (m.end()+m.start())/2 ) )
                if ( m.start() >= len(snp_alt)  and cur_suggested_shift > cur_excluded_seq_strat ):
                    cur_excluded_seq_strat =  cur_suggested_shift
            
            #print "cur_excluded_seq_strat:" + str(cur_excluded_seq_strat)
            
            # in case the end of the chromosome was reached
            if donor_nt_add_right > len(right_seq):
                print("shifting donor to avoid the end of the chromosome")
                cur_excluded_seq_strat = cur_excluded_seq_strat + (donor_nt_add_right-len(right_seq))
            
            
            #print "cur_excluded_seq_strat:" + str(cur_excluded_seq_strat)
            
            
            if cur_excluded_seq_strat>0:
                #print "moving to the left...."
                donor_seq_offset = donor_seq_offset - cur_excluded_seq_strat
                donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
                left_seq = str(genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ])
                donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq))
                right_seq = str( genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)])
        
                
                #print "2 donor_seq_offset:  %d" % (donor_seq_offset)
                
                # after moving left no problem on the right side
                cur_excluded_seq_strat = 0
            

                cur_excluded_seq_end = 0
                for m in re.finditer('|'.join(excluded_seqs),left_seq + str(Seq(snp_alt,generic_dna))  ):
                    #print m.end() 
                    if ( m.end() < len(left_seq)) and m.end() > cur_excluded_seq_end:
                        cur_excluded_seq_end = int( np.floor( (m.end()+m.start())/2 ) )
                       
                #print "cur_excluded_seq_end:" + str(cur_excluded_seq_end)
   
            
            if (cur_excluded_seq_strat == 0 and cur_excluded_seq_end == 0 and
                len(left_seq) >= min_dist_cut_to_donor_edge and len(right_seq) >= min_dist_cut_to_donor_edge):
                # found ok sequence
                donor_nt_add_left = len(left_seq)
                donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq)) # len(right_seq)
            else:
                warnings.warn("Can NOT get donor with no excluded sequences for guide_id: %s, length left: %d, length right %d (%s,%s) (%d,%d)" % (guide_id, len(left_seq), len(right_seq), left_seq, right_seq, cur_excluded_seq_strat, cur_excluded_seq_end) )
                return (out_guide_donor_df, donor_id)

            #print "found donor without excluded"
            #print int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut + 0 )
            #print donor_nt_add_left
            
            #print "3 donor_seq_offset:  %d" % (donor_seq_offset)
            #donor_seq_offset = int( np.floor((donor_length-len_diff_alt_ref)/2) + shift_dis_from_cut + 0 ) - donor_nt_add_left
            
        else:
            # num of nt to add left and right 
            donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
            donor_nt_add_right = int( (donor_length-len(snp_alt)) - donor_nt_add_left )
        
        
        #print "++++++"
        #print "snp_chr_pos %d" % (snp_chr_pos)
        #print "len(snp_ref) %d" % (len(snp_ref))
        #print "donor_nt_add_right %d" % (donor_nt_add_right)
        
        #print "DEBUG coor %d,  %d,  %d" % ( (snp_chr_pos - donor_nt_add_left) , (snp_chr_pos + len(snp_ref) + donor_nt_add_right), len(genome_seq[guide_chrom].seq))
        #print "DEBUG %s , %s ,  %s" % (str (genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ]) ,
        #                               str(Seq(snp_alt,generic_dna)),
        #                               str(genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)]))
        
        
        # if the SNP is at the edge of the chromosome give warning and do not return a donor        
        if ( (snp_chr_pos - donor_nt_add_left)  <  0 or (snp_chr_pos + len(snp_ref) + donor_nt_add_right) > len(genome_seq[guide_chrom].seq) ):
            warnings.warn("Can NOT get donor for SNP because on the edge of the chromosome. guide_id: %s, length left: %d, length right %d (%s,%s), minimal cooridnate: %d, maximal coordinate: %d, chromosome length: %d" % \
                          (guide_id, len(left_seq), len(right_seq), left_seq, right_seq, (snp_chr_pos - donor_nt_add_left), (snp_chr_pos + len(snp_ref) + donor_nt_add_right) , len(genome_seq[guide_chrom].seq) ) )
            return (out_guide_donor_df, donor_id)
            

            
        donor_seq = genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ] + \
                    Seq(snp_alt,generic_dna)  + \
                    genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)]

        #print "donor_seq before revcomp:  %s" % (str(donor_seq))
        # donor orientation should match the guide
        if (guide_is_negative_strand):
                donor_seq = donor_seq.reverse_complement()
                #print "revcomp the donor seq"
        
        #print "donor_seq after revcomp:  %s" % (str(donor_seq))
        
        
        
        # appending to donor sequences matrix
        donor_id += 1 
        
        cur_donor_line = pd.DataFrame({
            'var_id' : pd.Series(var_id), 'guide_id' : pd.Series(guide_id),
            'donor_id' : pd.Series(guide_id + ':' + donor_mut_type + ':offset' + str(donor_seq_offset) + ':donorID' + str(donor_id) + ':EditPosInGuide' +  str(SNP_pos_in_guide)), 
            'donor_seq': pd.Series(str(donor_seq)), 
            'donor_seq_shift' : pd.Series(int(donor_seq_offset)), 
            'donor_mut_pos_in_guide' : pd.Series(int(SNP_pos_in_guide)), 
            'donor_info_str' : pd.Series("ALT2REF:" + str(snp_ref) + ">" + str(snp_alt)),
            'set_name' : pd.Series(str(set_name))   })
        out_guide_donor_df = out_guide_donor_df.append(cur_donor_line)
        
    elif  (donor_mut_name == 'twoREF2ALT'):
        
        
        #print("In design of " + donor_mut_name)
        # TODO implement shifting to avoid sequences
        
        
        second_var_ids_str = str(snps_df.loc[snps_df['var_id'] == var_id, 'second_var_ids'].values[0])
        if second_var_ids_str == 'nan':
            second_var_ids = []
        else:
            second_var_ids =  second_var_ids_str.split(',')
        
        
        #print('#' * 20)
        #print(var_id)
        #print(second_var_ids)
        
        for second_var_id in second_var_ids:
            
            #print(second_var_id)
            cur_second_snp_df = snps_df.loc[snps_df['var_id'] == second_var_id, :]
            #print(cur_second_snp_df)
            
            
            second_snp_ref = str(cur_second_snp_df['REF'].values[0])
            second_snp_alt = (str((cur_second_snp_df['ALT'].values[0]))).replace('[','').replace(']','') # assuming '[]' was removed from the string
            second_snp_chr_pos = int(cur_second_snp_df['SNP_chr_pos'].values[0]) - 1 # TODO this is one based? (VCF)
            
            second_snp_chrom = str(cur_second_snp_df['chrom'].values[0])
            
            #  To be more careful you can change  warnings.warn, continue to raise ValueError
            
            if (guide_chrom != second_snp_chrom):
                warnings.warn(("both snps do not have the same chr %s, %s - %s, %s" %  (var_id,   second_var_id,
                                                                                           guide_chrom, second_snp_chrom) ))
                continue
            if (snp_chr_pos ==  second_snp_chr_pos):
                warnings.warn(("both snps have the same coordinate %s, %s - %d, %d" %  (var_id,   second_var_id,
                                                                                           snp_chr_pos, second_snp_chr_pos) ))
                continue
            elif (snp_chr_pos <  second_snp_chr_pos):
                if (snp_chr_pos + len(snp_ref) -1 >=  second_snp_chr_pos):
                    warnings.warn(("snp ref overlap (1) %s, %s - %d, %d - %s, %s" %  (var_id,   second_var_id,
                                                                                     snp_chr_pos, second_snp_chr_pos,
                                                                                     snp_ref, second_snp_ref) ))
                    continue
                
                
                two_snps_seq = ( Seq(snp_alt,generic_dna) +
                                 genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (second_snp_chr_pos)] +
                                 Seq(second_snp_alt,generic_dna) )
                left_coor = snp_chr_pos
                right_coor = second_snp_chr_pos + len(second_snp_ref)
                
            else: # (snp_chr_pos >  second_snp_chr_pos):
                if (second_snp_chr_pos + len(second_snp_ref) -1 >= snp_chr_pos):
                    warnings.warn(("snp ref overlap (2) %s, %s - %d, %d - %s, %s" %  (var_id,   second_var_id,
                                                                                     snp_chr_pos, second_snp_chr_pos,
                                                                                     snp_ref, second_snp_ref) ))
                    continue
                
                two_snps_seq = ( Seq(second_snp_alt,generic_dna) +
                                 genome_seq[guide_chrom].seq[ (second_snp_chr_pos + len(second_snp_ref)) : (snp_chr_pos)] +
                                 Seq(snp_alt,generic_dna) )
                left_coor = second_snp_chr_pos
                right_coor = snp_chr_pos + len(snp_ref)
                
            
            # adding up/downstream seqeucnes
            cur_donor_len = len(two_snps_seq)
            
            # distance between snps too long
            if (cur_donor_len > donor_length):
                warnings.warn(("donor length too long (1) %s, %s - %d, %d - %s, %s - %s" %  (var_id,   second_var_id,
                                                                                     snp_chr_pos, second_snp_chr_pos,
                                                                                     snp_ref, second_snp_ref,
                                                                                     str(two_snps_seq)) ))
                continue
            
            start_left_coor = min(left_coor,  guide0_chr_pos - min_dist_cut_to_donor_edge)
            end_right_coor =  max(right_coor, guide0_chr_pos + min_dist_cut_to_donor_edge + 1)
            
            if ( (end_right_coor -  start_left_coor) > donor_length):
                warnings.warn(("donor length too long (2) %s, %s - %d, %d - %s, %s - %s - %d, %d" %  (var_id,   second_var_id,
                                                                                     snp_chr_pos, second_snp_chr_pos,
                                                                                     snp_ref, second_snp_ref,
                                                                                     str(two_snps_seq),
                                                                                     start_left_coor, end_right_coor) ))
                continue
            
            add_to_right = int((donor_length - (end_right_coor -  start_left_coor)) /2)
            add_to_left = int(donor_length - (end_right_coor -  start_left_coor)) - add_to_right
            
            
            start_left_coor = start_left_coor - add_to_left
            end_right_coor = end_right_coor + add_to_right
            
            
            donor_seq = (genome_seq[guide_chrom].seq[ start_left_coor:left_coor] + 
                        two_snps_seq  + 
                        genome_seq[guide_chrom].seq[right_coor:end_right_coor])

            # donor orientation should match the guide
            if (guide_is_negative_strand):
                    donor_seq = donor_seq.reverse_complement()
            
            donor_seq_offset = 0
            
            dist_between_snps = int(np.abs(snp_chr_pos -  second_snp_chr_pos))
            
            donor_id += 1 
        
            cur_donor_line = pd.DataFrame({
                'var_id' : pd.Series(var_id),
                'guide_id' : pd.Series(guide_id),
                'donor_id' : pd.Series( ( guide_id + ':' + donor_mut_type +
                                         ':offset' + str(donor_seq_offset) +
                                         ':donorID' + str(donor_id) +
                                         ':EditPosInGuide' +  str(SNP_pos_in_guide) +
                                         ':FirstSNPID^' + str(var_id)  +
                                         ':SecondSNPID^' + str(second_var_id) +
                                         ':DistBetweenSNPs^' + str(dist_between_snps)
                                        )), 
                'donor_seq': pd.Series(str(donor_seq)), 
                'donor_seq_shift' : pd.Series(int(donor_seq_offset)), 
                'donor_mut_pos_in_guide' : pd.Series(int(SNP_pos_in_guide)), 
                'donor_info_str' : pd.Series("REF2ALT_twoSNPs:" + str(snp_ref) + ">" + str(snp_alt) + ":" + str(second_snp_ref) + ">" + str(second_snp_alt)),
                'set_name' : pd.Series(str(set_name))   })
            out_guide_donor_df = out_guide_donor_df.append(cur_donor_line)
            
        #print "returning donor"
    elif  (donor_mut_name == 'multiREF2ALT'):
        
        
        #print("In design of " + donor_mut_name)
        # TODO implement shifting to avoid sequences
        
       
        
        
        
        cur_snps_df = snps_df.loc[ snps_df['var_id'].isin(to_modify_var_ids),:].copy()
        
        cur_snps_df.sort_values(by = 'SNP_chr_pos',inplace = True, ascending = True)
        
        
        # makesure all are on the same chromosome 
        assert( np.all((cur_snps_df.chrom == cur_snps_df.chrom.values[0]).values))
        
        if (guide_chrom != cur_snps_df.chrom.values[0]):
                raise Value(("snps and guide do not have the same chr %s, %s - %s, %s" %  (var_id,   guide_id,
                                                                                           guide_chrom, cur_snps_df.chrom.values[0]) ))
        
    
        
        first_snp_chr_pos = cur_snps_df.SNP_chr_pos.values[0] - 1
        last_snp_chr_pos = cur_snps_df.SNP_chr_pos.values[-1] - 1
        
        
        edited_seq = Seq('',generic_dna)
        previous_snp_ref_end = None
        
        left_coor = first_snp_chr_pos
        right_coor = last_snp_chr_pos
        
        var_ids_str =''
        dist_between_snps = ''
        snps_ref2alt_str = ''
        
        for edited_var_id in cur_snps_df.var_id.values:

            #print(second_var_id)
            edited_snp_df = snps_df.loc[snps_df['var_id'] == edited_var_id, :]
            #print(cur_second_snp_df)
            
            
            # getting cur edited SNP details
            cur_snp_ref = str(edited_snp_df['REF'].values[0])
            cur_snp_alt = (str((edited_snp_df['ALT'].values[0]))).replace('[','').replace(']','') 
            cur_snp_chr_pos = int(edited_snp_df['SNP_chr_pos'].values[0]) - 1 # TODO this is one based? (VCF)
            cur_snp_chrom = str(edited_snp_df['chrom'].values[0])
            
            #  To be more careful you can change  warnings.warn, continue to raise ValueError
            #if (guide_chrom != cur_snp_chrom):
            #    raise Value(("both snps do not have the same chr %s, %s - %s, %s" %  (var_id,   second_var_id,
            #                                                                               guide_chrom, second_snp_chrom) ))
                
            
            if (previous_snp_ref_end is not None):
                edited_seq = (edited_seq +
                              genome_seq[guide_chrom].seq[ previous_snp_ref_end : cur_snp_chr_pos] +
                              Seq(cur_snp_alt,generic_dna)
                             )
                dist_between_snps = dist_between_snps + '^' +  str(int(cur_snp_chr_pos - previous_snp_ref_end + 1))
            else:
                edited_seq = Seq(cur_snp_alt,generic_dna)
                
            previous_snp_ref_end =  (cur_snp_chr_pos + len(cur_snp_ref))
            
            left_coor = min(left_coor, cur_snp_chr_pos)
            right_coor = max(right_coor, previous_snp_ref_end)
            
            var_ids_str = var_ids_str + '^' + edited_var_id
            snps_ref2alt_str = snps_ref2alt_str + cur_snp_ref + ">" + cur_snp_alt + "^"
            
            
        # adding up/downstream seqeucnes
        cur_donor_len = len(edited_seq)
        

        start_left_coor = min(left_coor  - min_dist_cut_to_donor_edge,     guide0_chr_pos - min_dist_cut_to_donor_edge)
        end_right_coor =  max(right_coor - min_dist_cut_to_donor_edge + 1, guide0_chr_pos + min_dist_cut_to_donor_edge + 1)
        
        
        donor_seq = (genome_seq[guide_chrom].seq[ start_left_coor:left_coor] + 
                    edited_seq  + 
                    genome_seq[guide_chrom].seq[right_coor:end_right_coor])
        

        if ( len(donor_seq) < donor_length):
            add_to_right = int( (donor_length - len(donor_seq)) / 2 )
            add_to_left = int(   donor_length - len(donor_seq)) - add_to_right
        
        
            start_left_coor_addition = start_left_coor - add_to_left
            end_right_coor_addition = end_right_coor + add_to_right
        
            donor_seq = (genome_seq[guide_chrom].seq[ start_left_coor_addition:start_left_coor] + 
                        donor_seq  + 
                        genome_seq[guide_chrom].seq[end_right_coor:end_right_coor_addition])

        # donor orientation should match the guide
        if (guide_is_negative_strand):
                donor_seq = donor_seq.reverse_complement()
        
        donor_seq_offset = 0
        
        if ( len(donor_seq) > donor_length):
            warnings.warn(("donor length too long (2) %s, %s - %d, %d %d - %s - %d, %d" %  (var_id,   guide_id,
                                                                                 left_coor, right_coor,guide0_chr_pos, 
                                                                                 str(edited_seq),
                                                                                 start_left_coor, end_right_coor) ))
        else:
            donor_id += 1 
        
            cur_donor_line = pd.DataFrame({
                'var_id' : pd.Series(var_id),
                'guide_id' : pd.Series(guide_id),
                'donor_id' : pd.Series( ( guide_id + ':' + donor_mut_type +
                                         ':offset' + str(donor_seq_offset) +
                                         ':donorID' + str(donor_id) +
                                         ':EditPosInGuide' +  str(SNP_pos_in_guide) +
                                         ':SNPIDs' + str(var_ids_str)  +
                                         ':DistBetweenSNPs' + str(dist_between_snps)
                                        )), 
                'donor_seq': pd.Series(str(donor_seq)), 
                'donor_seq_shift' : pd.Series(int(donor_seq_offset)), 
                'donor_mut_pos_in_guide' : pd.Series(int(SNP_pos_in_guide)), 
                'donor_info_str' : pd.Series("REF2ALT_multiSNPs:" + str(snps_ref2alt_str)),
                'set_name' : pd.Series(str(set_name))   })
            out_guide_donor_df = out_guide_donor_df.append(cur_donor_line)
            
        #print "returning donor"
        
    else:
        raise ValueError('get_donor_mut_for_SNP_guide unknown donor_mut_type:' + donor_mut_type +  " and donor_mut_name:" + donor_mut_name)
          

    return (out_guide_donor_df, donor_id)

###############################################################################################################
def get_donor_mut_for_SNP_guide(guide_id, snps_df, guides_df, genome_seq,
                                donor_mut_type, donor_length, excluded_seqs,
                                donor_seq_offset, min_dist_5prime_arm, min_dist_3prime_arm,
                                set_name, donor_id,
                                to_modify_var_ids = None):
    """
    simplified version of get_donor_mut_for_SNP_guide, only does REF2ALT donor design
    adjusted donor offset to account for guide strandedness (i.e. assymetric donor arms follow forward and reverse orientations)
    TO-DO: CONSIDER STRATEGY TO GIVE SUFFICIENT 3' HOMOLOGY FOR REVERSE GUIDES DOING INSERTIONS (e.g. fix short arm with lower limit first?)
    """
    
    out_guide_donor_df = pd.DataFrame(data=None)
    
    # select guide
    guides_df = guides_df[guides_df['guide_id'] == guide_id]

    # general
    guide_cut_chr_pos = int(guides_df['guide_cut_chr_pos'].iloc[0])
    guide0_chr_pos = int(guides_df['guide0_chr_pos'].iloc[0])
    guide_is_negative_strand = (guides_df['guide_strand'].iloc[0] == '-')
    guide_chrom = str(guides_df['chrom'].iloc[0])
    
    # specific to SNPs
    var_id = str(guides_df['var_id'].iloc[0])
    SNP_pos_in_guide = int(guides_df['SNP_pos_in_guide'].iloc[0])
    snp_ref = str(guides_df['REF'].iloc[0])
    snp_alt = str((guides_df['ALT'].iloc[0])) # assuming '[]' was removed from the string
    snp_chr_pos = int(guides_df['SNP_chr_pos'].iloc[0]) - 1 # gives 0-indexed position for SNP
    
    donor_mut_type_splt = donor_mut_type.split('_')
    donor_mut_name = donor_mut_type_splt[0]

    donor_seq = Seq("", generic_dna)
    
    
    if  (donor_mut_name == 'REF2ALT'):
        # shift donor by distance to Cas9 cut
        # positive means donor shifted left. negative means donor shifted right
        shift_dis_from_cut =  int(snp_chr_pos - guide_cut_chr_pos + 1) 
        
        # determine initial donor offset
        # if donor offset is not specified, donor will be symmetric around Cas9 cut site
        # if donor offset is positive, donor will have shorter left arm and longer right arm
        # if donor offset is negative, donor will have longer left arm and shorter right arm
        if donor_seq_offset is None:
            donor_seq_offset = 0
        # if guide is on negative strand, flip donor offset sign so that 5' and 3' offsets follow the same direction
        elif guide_is_negative_strand:
            donor_seq_offset = donor_seq_offset *-1
        
        
        # inspect donor sequence for excluded sequences
        # search donor left arm, shift sequence right till excluded sequences are removed
        donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
        left_seq = str(genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ]) 

        cur_excluded_seq_end = 0
        for m in re.finditer('|'.join(excluded_seqs),left_seq + str(Seq(snp_alt,generic_dna))  ):
            if ( m.end() < len(left_seq)) and m.end() > cur_excluded_seq_end:
                cur_excluded_seq_end = int( np.floor( (m.end()+m.start())/2 ) ) # half of rightmost excluded sequence
                    
        if cur_excluded_seq_end>0:
            donor_seq_offset = donor_seq_offset + cur_excluded_seq_end
            donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
            left_seq = str(genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ])
            # excluded sequences removed from left arm
            cur_excluded_seq_end = 0
        
        
        # search donor right arm, shift sequence left till excluded sequences are removed
        donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq))
        right_seq = str( genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)])

        cur_excluded_seq_start = 0
        for m in re.finditer('|'.join(excluded_seqs),str(Seq(snp_alt,generic_dna)) + right_seq):
            cur_suggested_shift = len(right_seq) - int( np.floor( (m.end()+m.start())/2 ) )
            if ( m.start() >= len(snp_alt)  and cur_suggested_shift > cur_excluded_seq_start ):
                cur_excluded_seq_start =  cur_suggested_shift

        # in case the end of the chromosome was reached
        if donor_nt_add_right > len(right_seq):
            print("shifting donor to avoid the end of the chromosome")
            cur_excluded_seq_start = cur_excluded_seq_start + (donor_nt_add_right-len(right_seq))

        if cur_excluded_seq_start>0:
            donor_seq_offset = donor_seq_offset - cur_excluded_seq_start
            donor_nt_add_left = int( np.floor((donor_length-len(snp_alt))/2) + shift_dis_from_cut - donor_seq_offset )
            left_seq = str(genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ])
            donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq))
            right_seq = str( genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)])
            # excluded sequences removed from right arm
            cur_excluded_seq_start = 0
            
            # inspect left arm again to see if excluded sequences could not be removed after adjusting right arm
            cur_excluded_seq_end = 0
            for m in re.finditer('|'.join(excluded_seqs),left_seq + str(Seq(snp_alt,generic_dna))  ):
                if ( m.end() < len(left_seq)) and m.end() > cur_excluded_seq_end:
                    cur_excluded_seq_end = int( np.floor( (m.end()+m.start())/2 ) )


        # calculate cut to donor edge arm lengths
        # distance calculations account for guide strand
        if (guide_is_negative_strand):
            dist_3prime_arm = int( len(left_seq) - shift_dis_from_cut + 1 ) 
            dist_5prime_arm = int( donor_length - dist_3prime_arm )
        else:
            dist_5prime_arm = int( len(left_seq) - shift_dis_from_cut + 1 )
            dist_3prime_arm = int( donor_length - dist_5prime_arm )
        
#         print(donor_id, guide_is_negative_strand, dist_5prime_arm, dist_3prime_arm)
        
        
        # do not return a donor if excluded sequences found
        if (cur_excluded_seq_start > 0) or (cur_excluded_seq_end > 0):
            warnings.warn("No valid donor for guide id %s: excluded sequence found (%d,%d). (left seq: %s, right seq: %s)" % (guide_id, cur_excluded_seq_start, cur_excluded_seq_end, left_seq, right_seq) )
            return (out_guide_donor_df, donor_id)
        
        # do not return a donor if donor arms do not meet minimum distance cutoffs
        if (dist_5prime_arm < min_dist_5prime_arm) or (dist_3prime_arm < min_dist_3prime_arm):
            warnings.warn("No valid donor for guide id %s: donor arm too short (5' arm length: %d (min req: %d), 3' arm length: %d (min req: %d)). (left seq: %s, right seq: %s)" % (guide_id, dist_5prime_arm, min_dist_5prime_arm, dist_3prime_arm, min_dist_3prime_arm, left_seq, right_seq) )
            return (out_guide_donor_df, donor_id)
        
        # if donor passes, donor sequence OK
        donor_nt_add_left = len(left_seq)
        donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq)) # len(right_seq)
        
#         # get finalized left and right donor arms if no excluded sequences and meets minimum donor arm distances
#         if (cur_excluded_seq_start == 0 and cur_excluded_seq_end == 0 and
#             dist_5prime_arm >= min_dist_5prime_arm and dist_3prime_arm >= min_dist_3prime_arm):
#             donor_nt_add_left = len(left_seq)
#             donor_nt_add_right = int( (donor_length-len(snp_alt)) - len(left_seq)) # len(right_seq)
#         # do not return a donor if excluded sequences found or distance to donor edge is below cutoff
#         else:
#             warnings.warn("Can NOT get donor with no excluded sequences for guide_id: %s, length left: %d, length right %d (%s,%s) (%d,%d)" % (guide_id, len(left_seq), len(right_seq), left_seq, right_seq, cur_excluded_seq_start, cur_excluded_seq_end) )
#             return (out_guide_donor_df, donor_id)

        
        # do not return a donor if the SNP is at chromosome edge    
        if ( (snp_chr_pos - donor_nt_add_left)  <  0 or (snp_chr_pos + len(snp_ref) + donor_nt_add_right) > len(genome_seq[guide_chrom].seq) ):
            warnings.warn("Can NOT get donor for SNP because on the edge of the chromosome. guide_id: %s, length left: %d, length right %d (%s,%s), minimal cooridnate: %d, maximal coordinate: %d, chromosome length: %d" % \
                          (guide_id, len(left_seq), len(right_seq), left_seq, right_seq, (snp_chr_pos - donor_nt_add_left), (snp_chr_pos + len(snp_ref) + donor_nt_add_right) , len(genome_seq[guide_chrom].seq) ) )
            return (out_guide_donor_df, donor_id)
        
        
        # assemble finalized donor sequence
        donor_seq = genome_seq[guide_chrom].seq[ (snp_chr_pos - donor_nt_add_left) : (snp_chr_pos) ] + \
                    Seq(snp_alt,generic_dna)  + \
                    genome_seq[guide_chrom].seq[ (snp_chr_pos + len(snp_ref)) : (snp_chr_pos + len(snp_ref) + donor_nt_add_right)]
        # donor orientation should match guide
        if (guide_is_negative_strand):
            donor_seq = donor_seq.reverse_complement()
        # appending to donor sequences matrix
        donor_id += 1 
        
        cur_donor_line = pd.DataFrame({
            'var_id' : pd.Series(var_id), 'guide_id' : pd.Series(guide_id),
            'donor_id' : pd.Series(guide_id + ':' + donor_mut_type + ':offset' + str(donor_seq_offset) + ':donorID' + str(donor_id) + ':EditPosInGuide' +  str(SNP_pos_in_guide)), 
            'donor_seq': pd.Series(str(donor_seq)), 
            'donor_seq_shift' : pd.Series(int(donor_seq_offset)), 
            'donor_mut_pos_in_guide' : pd.Series(int(SNP_pos_in_guide)), 
            'donor_info_str' : pd.Series("REF2ALT:" + str(snp_ref) + ">" + str(snp_alt)),
            'set_name' : pd.Series(str(set_name))   })
        out_guide_donor_df = out_guide_donor_df.append(cur_donor_line)

    return (out_guide_donor_df, donor_id)

    
###############################################################################################################
def design_donor_for_SNP_guides_old2020(input_design_table_filename, 
                                input_SNP_table_filename, input_guides_with_features_table_filename,
                                input_genome_fasta_filename,
                                donor_length, excluded_seqs, min_dist_cut_to_donor_edge,
                                output_donor_table_filename):
    
    ################################################################################################
    # loading input data frames
    ################################################################################################
    
    for excluded_seq in excluded_seqs:
        print("Excluded seq: %s" % (excluded_seq))
    
    
    # load snps df
    snps_df = pd.read_table(input_SNP_table_filename, sep='\t', na_values = "")

    # load guide df
    guides_df = pd.read_table(input_guides_with_features_table_filename, sep='\t', na_values = "")
    
    guides_df['ALT'] = guides_df['ALT'].str.replace('[','').str.replace(']','')
    
    #guides_df['ALT'] = guides_df['ALT'].apply(ast.literal_eval)

    # loading genome fasta file
    genome_seq = SeqIO.to_dict(SeqIO.parse(open(input_genome_fasta_filename),'fasta', alphabet=generic_dna))
    
    # loading design table
    donors_design_df = read_snp_donor_design_table(input_design_table_filename)


    out_guide_donor_df = pd.DataFrame(data=None)

    # running id to create unique id
    donor_id = -1 
    # iterating of the design matrix
    for idx,row in donors_design_df.iterrows():
        

        cur_guide_ids = get_design_filter_guide_ids(row, snps_df, guides_df)
        cur_mut_type = str(row['donor_mut_type']).strip()
        cur_donor_seq_offsets = row['donor_seq_offsets']
        cur_set_name = row['set_name']
        
        if (idx % 1 == 0):
            print("--------- desgining donors according to line: %d  (# guide ids = %d)----------" % (idx, len(cur_guide_ids)))
            print(row)
            print("-----------------------------------------------------------")

        cur_gene_donor_df, donor_id = get_donor_mut_for_SNP_guides(guide_ids = cur_guide_ids, 
                                                         snps_df = snps_df, guides_df = guides_df, 
                                                         genome_seq = genome_seq,
                                                         donor_mut_type = cur_mut_type, 
                                                         donor_length = donor_length,
                                                         excluded_seqs = excluded_seqs,
                                                         donor_seq_offsets = cur_donor_seq_offsets,
                                                         min_dist_cut_to_donor_edge = min_dist_cut_to_donor_edge,
                                                         set_name = cur_set_name,
                                                         donor_id = donor_id)

        #print cur_gene_donor_df.shape[0]

        out_guide_donor_df = pd.concat([out_guide_donor_df,cur_gene_donor_df],ignore_index=True)
        
        #print out_guide_donor_df.shape[0]


    out_guide_donor_df['contains_excluded_seqs'] = False
    if 'donor_seq' in out_guide_donor_df.columns:
        out_guide_donor_df['contains_excluded_seqs'] = out_guide_donor_df['donor_seq'].str.contains( '|'.join(excluded_seqs) )
    
    
    
    if len(output_donor_table_filename) > 0:
            print("saving donor sequences to: " +  output_donor_table_filename)
            out_guide_donor_df.to_csv(output_donor_table_filename, sep='\t', index = False)

    return(out_guide_donor_df)
    

    
###############################################################################################################
def design_donor_for_SNP_guides(input_design_table_filename, 
                                input_SNP_table_filename, input_guides_with_features_table_filename,
                                input_genome_fasta_filename,
                                donor_length, excluded_seqs, min_dist_5prime_arm, min_dist_3prime_arm,
                                output_donor_table_filename):
    """
    modified to accept asymmetrical minimum cut to donor edge distances
    """

    ################################################################################################
    # loading input data frames
    ################################################################################################
    
    for excluded_seq in excluded_seqs:
        print("Excluded seq: %s" % (excluded_seq))
    
    
    # load snps df
    snps_df = pd.read_table(input_SNP_table_filename, sep='\t', na_values = "")

    # load guide df
    guides_df = pd.read_table(input_guides_with_features_table_filename, sep='\t', na_values = "")
    
    guides_df['ALT'] = guides_df['ALT'].str.replace('[','').str.replace(']','')
    
    #guides_df['ALT'] = guides_df['ALT'].apply(ast.literal_eval)

    # loading genome fasta file
    genome_seq = SeqIO.to_dict(SeqIO.parse(open(input_genome_fasta_filename),'fasta', alphabet=generic_dna))
    # set genome to uppercase
    genome_seq = {chrom : sequence.upper() for chrom, sequence in genome_seq.items()}
    
    # loading design table
    donors_design_df = read_snp_donor_design_table(input_design_table_filename)


    out_guide_donor_df = pd.DataFrame(data=None)

    # running id to create unique id
    donor_id = 0 
    # iterating of the design matrix
    for idx,row in donors_design_df.iterrows():
        

        cur_guide_ids = get_design_filter_guide_ids(row, snps_df, guides_df)
        cur_mut_type = str(row['donor_mut_type']).strip()
        cur_donor_seq_offsets = row['donor_seq_offsets']
        cur_set_name = row['set_name']
        
        if (idx % 1 == 0):
            print("--------- designing donors according to line: %d  (# guide ids = %d)----------" % (idx, len(cur_guide_ids)))
            print(row)
            print("-----------------------------------------------------------")

        cur_gene_donor_df, donor_id = get_donor_mut_for_SNP_guides(guide_ids = cur_guide_ids, 
                                                         snps_df = snps_df, guides_df = guides_df, 
                                                         genome_seq = genome_seq,
                                                         donor_mut_type = cur_mut_type, 
                                                         donor_length = donor_length,
                                                         excluded_seqs = excluded_seqs,
                                                         donor_seq_offsets = cur_donor_seq_offsets,
                                                         min_dist_5prime_arm = min_dist_5prime_arm, 
                                                         min_dist_3prime_arm = min_dist_3prime_arm,
                                                         set_name = cur_set_name,
                                                         donor_id = donor_id)

        #print cur_gene_donor_df.shape[0]

        out_guide_donor_df = pd.concat([out_guide_donor_df,cur_gene_donor_df],ignore_index=True)
        
        #print out_guide_donor_df.shape[0]


    out_guide_donor_df['contains_excluded_seqs'] = False
    if 'donor_seq' in out_guide_donor_df.columns:
        out_guide_donor_df['contains_excluded_seqs'] = out_guide_donor_df['donor_seq'].str.contains( '|'.join(excluded_seqs) )
    
    
    
    if len(output_donor_table_filename) > 0:
            print("saving donor sequences to: " +  output_donor_table_filename)
            out_guide_donor_df.to_csv(output_donor_table_filename, sep='\t', index = False)

    return(out_guide_donor_df)

###############################################################################################################
def design_donor_for_multipleSNP(input_multiSNP_design_table_filename, 
                                input_SNP_table_filename,
                                input_guides_with_features_table_filename,
                                input_genome_fasta_filename,
                                donor_length, excluded_seqs, min_dist_cut_to_donor_edge,
                                output_donor_table_filename):
    """
    input_guides_with_features_table_filename must be guides with featuere and ranked 
    """
    
    ################################################################################################
    # loading input data frames
    ################################################################################################
    
    for excluded_seq in excluded_seqs:
        print("Excluded seq: %s" % (excluded_seq))
    
    
    # load snps df
    snps_df = pd.read_table(input_SNP_table_filename, sep='\t', na_values = "")

    # load guide df
    guides_df = pd.read_table(input_guides_with_features_table_filename, sep='\t', na_values = "")
    
    # removing square brackets if exists
    guides_df['ALT'] = guides_df['ALT'].str.replace('[','').str.replace(']','')
    
    #guides_df['ALT'] = guides_df['ALT'].apply(ast.literal_eval)

    # loading genome fasta file
    genome_seq = SeqIO.to_dict(SeqIO.parse(open(input_genome_fasta_filename),'fasta', alphabet=generic_dna))
    
    # loading design table
    donors_design_df = read_snp_donor_design_table(input_multiSNP_design_table_filename)


    out_donor_df = pd.DataFrame(data=None)

    # running id to create unique id
    donor_id = -1 
    # iterating of the design matrix
    for idx,row in donors_design_df.iterrows():
        
        
        cur_mut_type = str(row['donor_mut_type']).strip()
        cur_donor_seq_offsets = row['donor_seq_offsets']
        cur_set_name = row['set_name']
        
        # for single snp no need to add each row
        # for pair of snps the snps table should contain 'second_var_ids' column.
        
        if (cur_mut_type in ['twoREF2ALT', 'REF2ALT']):
            
            cur_guide_ids = get_design_filter_guide_ids(row, snps_df, guides_df)

            
            if (idx % 1== 0):
                print("--------- desgining donors according to line: %d  (# guide ids = %d)----------" % (idx, len(cur_guide_ids)))
                print(row)
                print("-----------------------------------------------------------")
    
            cur_donors_df, donor_id = get_donor_mut_for_SNP_guides(guide_ids = cur_guide_ids, 
                                                             snps_df = snps_df, guides_df = guides_df, 
                                                             genome_seq = genome_seq,
                                                             donor_mut_type = cur_mut_type, 
                                                             donor_length = donor_length,
                                                             excluded_seqs = excluded_seqs,
                                                             donor_seq_offsets = cur_donor_seq_offsets,
                                                             min_dist_cut_to_donor_edge = min_dist_cut_to_donor_edge,
                                                             set_name = cur_set_name,
                                                             donor_id = donor_id)
        elif (cur_mut_type == 'multiREF2ALT'):
            
            if (idx % 100 == 0):
                print("--------- desgining donors according to line: %d ----------" % (idx))
                print(row)
                print("-----------------------------------------------------------")
                
            cur_target_var_id = row['target_var_id']
            cur_max_num_of_guides = int(row['Max_number_of_guides'])
            cur_modify_var_ids = row['modify_var_ids'].split(',')
            
            # TODO - consider just adding it to the set
            if(cur_target_var_id not in cur_modify_var_ids):
                raise ValueError("target SNP is not in the edited SNPs list: %s, %s" % (cur_target_var_id, ','.join(cur_modify_var_ids)))
            
            
            # TODO - not implemented
            # donor_seq_offsets
            
            # get guides
            guide_ids = guides_df.loc[((guides_df['var_id'] == cur_target_var_id) &
                           (guides_df['quality_rank']<cur_max_num_of_guides)), 'guide_id'].values
            
            
            cur_donors_df =  pd.DataFrame(data=None)
    
            # design per guide
            for guide_id in guide_ids:
                
                # not relevant columns: filter_in filter_out
                cur_guide_donors_df, donor_id = get_donor_mut_for_SNP_guides(guide_id,
                        snps_df = snps_df, guides_df = guides_df, genome_seq = genome_seq,
                        donor_mut_type = cur_mut_type,
                        donor_length = donor_length,
                        excluded_seqs = excluded_seqs, donor_seq_offset = None,
                        min_dist_cut_to_donor_edge = min_dist_cut_to_donor_edge,
                        set_name = cur_set_name, donor_id = donor_id,
                        to_modify_var_ids = cur_modify_var_ids)
                
                cur_donors_df = cur_donors_df.append(cur_guide_donors_df,ignore_index=True)            
        
        else:
            raise ValueError("unknown cur_mut_type: %s" % (cur_mut_type))
        

        # concating the cur table 
        out_donor_df = pd.concat([out_donor_df,cur_donors_df],ignore_index=True)
        
        #print out_guide_donor_df.shape[0]


    out_donor_df['contains_excluded_seqs'] = False
    if 'donor_seq' in out_donor_df.columns:
        out_donor_df['contains_excluded_seqs'] = out_donor_df['donor_seq'].str.contains( '|'.join(excluded_seqs) )
    
    
    
    if len(output_donor_table_filename) > 3:
            print("saving donor sequences to: " +  output_donor_table_filename)
            out_donor_df.to_csv(output_donor_table_filename, sep='\t', index = False)

    return(out_donor_df)
 



###############################################################################################################
def prepare_VEP_input_from_SNPs_df(snps_df_filename, vep_input_filename):
    snps_df = pd.read_table(snps_df_filename, sep='\t', na_values = "")
    annotate_df = snps_df.loc[snps_df.chrom != 'chr_mito', ['oligo_seq', 'var_id', 'chrom', 'SNP_chr_pos', 'REF', 'ALT', 'donor_info_str'] ].copy()

    annotate_df['ALT'] = annotate_df['ALT'].str.replace('[','').str.replace(']','')
    annotate_df['coor_start'] = annotate_df['SNP_chr_pos']
    annotate_df['coor_end'] = annotate_df['SNP_chr_pos'] + annotate_df['REF'].apply(lambda s: len(s)) - 1
    if annotate_df.chrom.str.isalpha().all():
        annotate_df['chrom_roman'] = annotate_df.chrom.str.replace('chr','')
    else:
        annotate_df['chrom_roman'] = annotate_df.chrom.apply(yeast_chr_numeric2roman).str.replace('chr','')
    annotate_df['strand'] = '+'

    annotate_df.loc[annotate_df.REF == '','REF'] = '-'
    annotate_df.loc[annotate_df.ALT == '','ALT'] = '-'
    annotate_df['allele_edit'] = annotate_df['REF'] + '/' + annotate_df['ALT']


    vep_input_df = annotate_df.loc[:,
                                   ['chrom_roman','coor_start', 'coor_end',
                                    'allele_edit', 'strand', 'var_id']]

    vep_input_df = vep_input_df.loc[np.isfinite(vep_input_df.coor_end),:]
    vep_input_df.ix[:,'coor_start'] = vep_input_df.coor_start.astype(int)
    vep_input_df.ix[:,'coor_end'] = vep_input_df.coor_end.astype(int)
    
    # writing VCF file
    vep_input_df.to_csv(vep_input_filename,sep='\t',index=False, header=False)

    return(vep_input_df)


###############################################################################################################
def rank_guides(cur_snp_guide_df, 
                off_targets_min_mismatch_SNP_guides = 10, 
                min_ok_Azimuth_score_SNP_guides=0.2, 
                edit_max_distance_from_PAM5prime = 18):
    
    # calculating if guide pass filter
    cur_snp_guide_df.loc[:,'is_passed_filter'] = ((cur_snp_guide_df['is_found_donor_with_no_excluded_seqs'].astype('bool').values == True) & 
                                                  (cur_snp_guide_df['Azimuth'].astype('float').values >= min_ok_Azimuth_score_SNP_guides) & 
                                                  (cur_snp_guide_df['SNP_pos_in_guide'].astype('int').values >= -edit_max_distance_from_PAM5prime) & 
                                                  (cur_snp_guide_df['is_nuclear_chromosome'].astype('bool').values == True))
    

    if (off_targets_min_mismatch_SNP_guides > 0):
        cur_snp_guide_df['is_passed_filter'] = cur_snp_guide_df['is_passed_filter'] & (cur_snp_guide_df['guide_map_mismatch_0'] == 0)
    if (off_targets_min_mismatch_SNP_guides > 1):
        cur_snp_guide_df['is_passed_filter'] = cur_snp_guide_df['is_passed_filter'] & (cur_snp_guide_df['guide_map_mismatch_1'] == 0)
    if (off_targets_min_mismatch_SNP_guides > 2):
        cur_snp_guide_df['is_passed_filter'] = cur_snp_guide_df['is_passed_filter'] & (cur_snp_guide_df['guide_map_mismatch_2'] == 0)
    if (off_targets_min_mismatch_SNP_guides > 3):
        cur_snp_guide_df['is_passed_filter'] = cur_snp_guide_df['is_passed_filter'] & (cur_snp_guide_df['guide_map_mismatch_3'] == 0)
        
    
    # ranking the guides
    # by filter, position, azimuth score, off targets

    tmp_cur_snp_guide_df = cur_snp_guide_df[  ['is_passed_filter', 'SNP_pos_in_guide','Azimuth', 'guide_map_mismatch_0','guide_map_mismatch_1','guide_map_mismatch_2','guide_map_mismatch_3'] ].copy()
    tmp_cur_snp_guide_df['rank_ind'] = range(tmp_cur_snp_guide_df.shape[0])
    tmp_cur_snp_guide_df.sort_values(['is_passed_filter', 'SNP_pos_in_guide', 'guide_map_mismatch_0','guide_map_mismatch_1','guide_map_mismatch_2','guide_map_mismatch_3','Azimuth'], 
              ascending=[False,False, True,True,True,True, False], inplace=True, na_position='last')

    cur_snp_guide_df.iloc[tmp_cur_snp_guide_df['rank_ind'].astype(int).values,cur_snp_guide_df.columns.get_loc('quality_rank')] = range(cur_snp_guide_df.shape[0])

    return(cur_snp_guide_df)

###############################################################################################################
def rank_and_filter_SNP_guides(input_guides_with_features_table_filename,
                              input_SNP_table_filename,
                              input_donor_table_filename = None,
                              output_guides_with_features_and_rank_table_filename = '',
                              off_targets_min_mismatch_SNP_guides = 10, 
                              min_ok_Azimuth_score_SNP_guides = 0.2, 
                              edit_max_distance_from_PAM5prime = 18):
    """
    For guides that target SNPs ranks the guides for each SNP ('quality_rank') and add a column marking whether it passed a filter
    """

    snp_guides_df = pd.read_table(input_guides_with_features_table_filename, sep='\t', na_values = "None")
    snp_df = pd.read_table(input_SNP_table_filename, sep='\t', na_values = "None")
    
    
    # columns required for the filtering
    if (input_donor_table_filename is not None):
        snp_donors_df = pd.read_table(input_donor_table_filename, sep='\t', na_values = "None")
        snp_guides_df['is_found_donor_with_no_excluded_seqs'] = snp_guides_df['guide_id'].isin( snp_donors_df['guide_id'][snp_donors_df['contains_excluded_seqs']==False] )
    else:
        warnings.warn("rank_and_filter_SNP_guides: donor table is not provided so does not filter according to excluded sequences in the donors (all get is_found_donor_with_no_excluded_seqs = True)"  )
        snp_guides_df['is_found_donor_with_no_excluded_seqs'] = True

    snp_guides_df['no_perfect_offtargets'] =  (snp_guides_df['guide_map_mismatch_0'] == 0)

    snp_guides_df['no_offtargets'] =  ((snp_guides_df['guide_map_mismatch_0'] == 0) &
                                       (snp_guides_df['guide_map_mismatch_1'] == 0) &
                                       (snp_guides_df['guide_map_mismatch_2'] == 0) &
                                       (snp_guides_df['guide_map_mismatch_3'] == 0) )


    snp_guides_df['is_passed_filter'] = False
    snp_guides_df['quality_rank'] = 1000

    # ranking guides for each SNP and adding a filter column
    snp_guides_df = snp_guides_df.groupby('var_id').apply(rank_guides, 
                                      off_targets_min_mismatch_SNP_guides = off_targets_min_mismatch_SNP_guides, 
                                      min_ok_Azimuth_score_SNP_guides=min_ok_Azimuth_score_SNP_guides, 
                                      edit_max_distance_from_PAM5prime = edit_max_distance_from_PAM5prime)

    # saving the table with the ranking and filter column
    if output_guides_with_features_and_rank_table_filename:
        print("Saving : " + output_guides_with_features_and_rank_table_filename)
        snp_guides_df.to_csv(output_guides_with_features_and_rank_table_filename, sep='\t', index = False)
  
    return(snp_guides_df)



###############################################################################################################
# functions for Bloom QTL merging 
###############################################################################################################


###############################################################################################################
def merge_SNPs_with_QTLs(input_SNP_table_filename, input_QTLs_filename, input_QTLs_H_filename,
                        input_selected_QTLs, input_known_genes_with_large_effect,
                        output_selected_qtls_table_filename, output_snp_in_selected_QTLs_filename,
                        output_SNP_merged_with_QTL_table_filename, output_SNP_withQTL_table_filename):
    """
    Merging SNPs with Bloom QTLs -
    output_SNP_withQTL_table_filename will contain additional columns and a 'set_name' column that allows to select SNPs for forther analysis
    """
    SNP_withQTLs_df = pd.DataFrame(data=None)
    
    ###########################################
    # loading and parsing the input tables
    ###########################################
    
    snps_df = pd.read_table(input_SNP_table_filename, sep='\t', na_values = "")

    qtls_df = pd.read_csv(input_QTLs_filename)
    qtls_H_df = pd.read_table(input_QTLs_H_filename, sep='\t', na_values = "None")

    selected_qtls_df = pd.read_table(input_selected_QTLs, sep='\t')
    gene_with_known_large_effect_df = pd.read_table(input_known_genes_with_large_effect, sep='\t')
    
    # parse - this is specific to Bloom QTLs
    chr_int2str = lambda x: "chr" + str(x).zfill(2)

    qtls_df['chrom'] = qtls_df['Chromosome'].apply(chr_int2str)
    qtls_df['qtl_start_POS'] = qtls_df['1.5 LOD Confidence Interval (Left) (bp)']
    qtls_df['qtl_end_POS'] = qtls_df['1.5 LOD Confidence Interval (Right) (bp)']
 
    # printing
    
    print("Selected traits:")
    print(selected_qtls_df)
    print("Genes with known effects:")
    print(gene_with_known_large_effect_df)
    
    ###########################################
    # select qtls 
    ###########################################

    is_qtl_not_contain_gene_with_known_effect = lambda x: not any([gene in str(x) for gene in gene_with_known_large_effect_df['Gene_name']])

    # selecting QTLs in Bloom's 8 (151) and removing 30 QTLs that do not contain genes with large effects (121)
    idx = qtls_df['Trait'].isin(selected_qtls_df['Trait']) & \
          qtls_df['Genes Under Confidence Interval (Standard Name)'].apply(is_qtl_not_contain_gene_with_known_effect)

    # selected QTLs table 
    sel_qtls_df = qtls_df[idx]
    print("Selected %d QTLs out of %d QTLs (filtering for selected QTLS and no known gene with large effect)" % (sel_qtls_df.shape[0],qtls_df.shape[0]))

    
    # writing merge between SNPs and selected QTLs. Notice that a SNP can be in more than one QTL
    if output_selected_qtls_table_filename:
        print("Saving : " + output_selected_qtls_table_filename)
        sel_qtls_df.to_csv(output_selected_qtls_table_filename, sep='\t', index = False)
        
 
    # TODO saving the excluded QTLs - add an input argument and add to calling scripts
    if False:
        
        is_qtl_contain_gene_with_known_effect = lambda x: any([gene in str(x) for gene in gene_with_known_large_effect_df['Gene_name']])
        
        # excluding QTLs in Bloom's 8 (151) and removing 30 QTLs that contain genes with large effects (121)
        idx = qtls_df['Trait'].isin(selected_qtls_df['Trait']) & \
              qtls_df['Genes Under Confidence Interval (Standard Name)'].apply(is_qtl_contain_gene_with_known_effect)
    
        # excluded QTLs table 
        excluded_qtls_df = qtls_df[idx]
        print("Excluded %d QTLs out of %d QTLs (filtering for selected QTLS and no known gene with large effect)" % (excluded_qtls_df.shape[0],qtls_df.shape[0]))
    
        # writing merge between SNPs and selected QTLs. Notice that a SNP can be in more than one QTL
        if output_excluded_qtls_table_filename:
            print("Saving : " + output_excluded_qtls_table_filename)
            excluded_qtls_df.to_csv(output_excluded_qtls_table_filename, sep='\t', index = False)
          
    ###########################################
    # merge the selected QTLs with 
    ###########################################
    snps_merged_with_QTLs_df = pd.merge(snps_df,sel_qtls_df,on='chrom')
    snps_merged_with_QTLs_df = snps_merged_with_QTLs_df[ (snps_merged_with_QTLs_df['qtl_start_POS'] <= snps_merged_with_QTLs_df['SNP_chr_pos']) & 
                                                 (snps_merged_with_QTLs_df['SNP_chr_pos']   <= snps_merged_with_QTLs_df['qtl_end_POS']) ]

    
    if output_SNP_merged_with_QTL_table_filename:
        print("Saving : " + output_SNP_merged_with_QTL_table_filename)
        snps_merged_with_QTLs_df.to_csv(output_SNP_merged_with_QTL_table_filename, sep='\t', index = False)
  
    
    ###########################################
    # unique SNPs in selected QTLs
    ###########################################
    
    #var_ids_in_selected_QTLs = snps_merged_with_QTLs_df['var_id'].unique()
    #var_ids_in_selected_QTLs_df = pd.DataFrame({'var_id':  pd.Series(var_ids_in_selected_QTLs)})
    
    #DEBUG
    #print(snps_merged_with_QTLs_df.shape[0] >0)
    #print(snps_merged_with_QTLs_df.groupby('var_id')[[ 'var_id','Trait' ] ].apply(lambda x: "%s" % ','.join(x['Trait'])))
    
    if (snps_merged_with_QTLs_df.shape[0] > 0):
        snps_traits_df = snps_merged_with_QTLs_df.groupby('var_id')[[ 'var_id','Trait' ] ].apply(lambda x: "%s" % ','.join(x['Trait'])).to_frame("Traits")
        snps_traits_df['var_id'] = snps_traits_df.index

        snps_traits2_df = snps_merged_with_QTLs_df.groupby('var_id')[[ 'var_id','Trait' ] ].size().to_frame("NumTraits")
        snps_traits2_df['var_id'] = snps_traits2_df.index

        snp_in_selected_QTLs_df = pd.merge(snps_traits_df, snps_traits2_df)
    else:
        print("WARNING snp_in_selected_QTLs_df is empty")
        snp_in_selected_QTLs_df = pd.DataFrame(None)

    
    
    print("#SNPs in seleced QTLs: " + str(snp_in_selected_QTLs_df.shape[0]) + "=" + str(len(snps_merged_with_QTLs_df['var_id'].unique())) + "(out of a total of " + str(len(snps_df['var_id'].unique())) + " SNPs)")

        
    # unique SNPIDs in selected QTLs
    if output_snp_in_selected_QTLs_filename:
        print("Saving : " + output_snp_in_selected_QTLs_filename)
        snp_in_selected_QTLs_df.to_csv(output_snp_in_selected_QTLs_filename, sep='\t', index = False)

    
    #print snp_in_selected_QTLs_df.head(5)
    
    ###########################################
    # adding QTL columns to the SNPs table
    ###########################################
    
    snps_withQTLs_df = snps_df.copy()
    if (snp_in_selected_QTLs_df.shape[0]>0):
        snps_withQTLs_df['is_in_selected_QTL'] = snps_withQTLs_df['var_id'].isin(snp_in_selected_QTLs_df['var_id'])
        snps_withQTLs_df = pd.merge(snps_withQTLs_df, snp_in_selected_QTLs_df, how='left', on=['var_id'])
    else:
        snps_withQTLs_df['is_in_selected_QTL'] = False
        snps_withQTLs_df['Traits'] = ''
        snps_withQTLs_df['NumTraits'] = 0                 
    
    snps_withQTLs_df['set_name'] = "not_in_sel_qtls"
    snps_withQTLs_df.ix[snps_withQTLs_df['is_in_selected_QTL'], 'set_name'] = 'selected_qtls'
    
    if output_SNP_withQTL_table_filename:
        print("Saving : " + output_SNP_withQTL_table_filename)
        snps_withQTLs_df.to_csv(output_SNP_withQTL_table_filename, sep='\t', index = False)
  

    return(snps_withQTLs_df)


###############################################################################################################
def generate_scanning_pseudo_snps(qtl_effect_df, row_index):
    """
    """
    cur_chr = qtl_effect_df['CHROM'][row_index]
    cur_start = qtl_effect_df['start_POS'][row_index]
    cur_end   = qtl_effect_df['end_POS'][row_index]
    
    chr_poistion_sr = pd.Series( np.arange(cur_start,cur_end+1))

    cur_df = pd.DataFrame({'POS' : chr_poistion_sr})
    cur_df['CHROM'] = cur_chr
    cur_df['Trait'] = qtl_effect_df['Trait'][row_index]
    cur_df['start_POS'] = qtl_effect_df['start_POS'][row_index]
    cur_df['end_POS'] = qtl_effect_df['end_POS'][row_index]
    
    return(cur_df)


###############################################################################################################
def generate_qtl_gene_guide_design_table(input_selected_qtls_table_filename, input_qtl_gene_design_table_filename,input_annot_gff_filename, include_dubious=False,
                                   output_qtl_uniq_gene_guides_design_table_filename = '', output_qtl_nonuniq_gene_guides_design_table_filename = '',
                                   output_qtl_genes_that_were_filtered_out_filename = ''):
    
    print("Generating design tables of guides for genes in selected QTLs")
    
    
    sel_qtls_df = pd.read_table(input_selected_qtls_table_filename, sep='\t', na_values = "")
    input_design_df = pd.read_table(input_qtl_gene_design_table_filename, sep='\t', na_values = "")
    input_design_df['key'] = 1 # for cartesian product
    
    
    genes_gff_df = sgd_gff2dataframe(input_annot_gff_filename, ['CDS'],include_dubious=include_dubious)
    
    
    nonuniq_qtl_genes_design_df = pd.DataFrame(data=None)

    for idx,row in sel_qtls_df.iterrows():
        cur_gene_ids = str(row['Genes Under Confidence Interval (Standard Name)']).split('|')
        #print cur_gene_ids
        cur_out_qtl_df = pd.DataFrame({"gene_id" : cur_gene_ids})
        cur_out_qtl_df['set_name'] = row['Trait']
        cur_out_qtl_df['Trait'] = row['Trait']
        cur_out_qtl_df['chrom'] = row['chrom']
        cur_out_qtl_df['qtl_start_POS'] = row['qtl_start_POS']
        cur_out_qtl_df['qtl_end_POS'] = row['qtl_end_POS']
        cur_out_qtl_df['Fraction of Phenotypic Variance Explained'] = row['Fraction of Phenotypic Variance Explained']
        cur_out_qtl_df['LOD'] = row['LOD']
        cur_out_qtl_df['QTL'] = row['Trait'] + "_" + row['chrom'] + str(row['qtl_start_POS']) + "_" + str(row['qtl_end_POS'])

        nonuniq_qtl_genes_design_df = pd.concat([nonuniq_qtl_genes_design_df,cur_out_qtl_df],ignore_index=True)


    print("Number of genes in QTLs (non unique, before filtering): %d" % (nonuniq_qtl_genes_design_df.shape[0]))


    # TODO -  add a filter criteria?
    filt_qtl_genes_design_df = nonuniq_qtl_genes_design_df.copy()
    
    # filter genes that do not appear in the gff

    filt_qtl_genes_design_df = filt_qtl_genes_design_df[filt_qtl_genes_design_df['gene_id'].isin(genes_gff_df['Name'])]
    
    
    if output_qtl_genes_that_were_filtered_out_filename:
        print("Saving : " + output_qtl_genes_that_were_filtered_out_filename)
        pd.DataFrame({"gene_id" : nonuniq_qtl_genes_design_df['gene_id'][~ nonuniq_qtl_genes_design_df['gene_id'].isin(filt_qtl_genes_design_df['gene_id'])].unique() }).to_csv(output_qtl_genes_that_were_filtered_out_filename, sep='\t', index = False)
        
    
    
    
    print("Number of genes in QTLs (non unique, filttered): %d" % (filt_qtl_genes_design_df.shape[0]))


    uniq_qtl_genes_design_df = filt_qtl_genes_design_df.groupby('gene_id') \
        [['Trait','gene_id'] ].apply(lambda x: "%s" % ','.join(x['Trait'])).to_frame("Traits")
    uniq_qtl_genes_design_df.reset_index(inplace=True)
    uniq_qtl_genes_design_df['set_name'] = 'selected_QTL_genes'
    
    
    print("Number of genes in QTLs (unique and filtered): %d" % (uniq_qtl_genes_design_df.shape[0]))

    # joining with design table
    nonuniq_qtl_genes_design_df['key'] = 1
    nonuniq_qtl_genes_design_df = pd.merge(nonuniq_qtl_genes_design_df, input_design_df, on='key')
    nonuniq_qtl_genes_design_df.drop('key', axis=1, inplace=True)


    uniq_qtl_genes_design_df['key'] = 1
    uniq_qtl_genes_design_df = pd.merge(uniq_qtl_genes_design_df, input_design_df, on='key')
    uniq_qtl_genes_design_df.drop('key', axis=1, inplace=True)
    
    if output_qtl_uniq_gene_guides_design_table_filename:
        print("Saving : " + output_qtl_uniq_gene_guides_design_table_filename)
        uniq_qtl_genes_design_df.to_csv(output_qtl_uniq_gene_guides_design_table_filename, sep='\t', index = False)
    
    if output_qtl_nonuniq_gene_guides_design_table_filename:
        print("Saving : " + output_qtl_nonuniq_gene_guides_design_table_filename)
        nonuniq_qtl_genes_design_df.to_csv(output_qtl_nonuniq_gene_guides_design_table_filename, sep='\t', index = False)


    return (uniq_qtl_genes_design_df, nonuniq_qtl_genes_design_df)

###############################################################################################################


###############################################################################################################
# functions for essential genes analysis 
###############################################################################################################

###############################################################################################################
def parse_essential_genes_to_design_table(essential_genes_filename,
                                        input_annot_gff_filename,
                                        input_selected_qtls_table_filename,
                                        qtl_trait_for_essential_genes,
                                        include_dubious,
                                        number_of_genes_to_select,
                                        shuffle_genes,
                                        input_design_table_filename,
                                        output_essential_gene_design_table_filename,
                                        shuffle_genes_rand_seed = 14):
    
    essential_genes_df = read_essential_genes(essential_genes_filename,
                                            input_annot_gff_filename,
                                            input_selected_qtls_table_filename,
                                            qtl_trait_for_essential_genes,
                                            include_dubious,
                                            number_of_genes_to_select,
                                            shuffle_genes,
                                            shuffle_genes_rand_seed = shuffle_genes_rand_seed)
    
    design_df = pd.read_table(input_design_table_filename, sep='\t', na_values = "None")
    design_df['donor_seq_offsets'] = design_df['donor_seq_offsets'].astype(str)
    design_df['donor_seq_offsets'].replace(to_replace="nan",value="[None]",inplace=True)
    
    essential_genes_df['key'] = 1
    design_df['key'] = 1
    
    num_essential_genes = essential_genes_df.shape[0]
    
    nonuniq_qtl_genes_design_df = pd.DataFrame(data=None)
    
    for idx,row in design_df.iterrows():
        
        #print "parsing raw %d" % (idx)
        
        if np.isnan(row['gene_ind_start']):
            gene_ind_start = 0    
        else:
            gene_ind_start = int(row['gene_ind_start'])
        
        if np.isnan(row['gene_ind_end']):
            gene_ind_end = num_essential_genes
        else:
            gene_ind_end = int(row['gene_ind_end'])
        
        cur_nonuniq_qtl_genes_design_df = pd.merge(essential_genes_df.iloc[gene_ind_start:gene_ind_end,:], (pd.DataFrame(row)).transpose(), on='key')
        
        nonuniq_qtl_genes_design_df =  pd.concat([nonuniq_qtl_genes_design_df,cur_nonuniq_qtl_genes_design_df],ignore_index=True)
        
    nonuniq_qtl_genes_design_df.drop('key', axis=1, inplace=True)
    
    
    nonuniq_qtl_genes_design_df.reset_index(inplace=True)
    
    if output_essential_gene_design_table_filename:
        print("Saving : " + output_essential_gene_design_table_filename)
        nonuniq_qtl_genes_design_df.to_csv(output_essential_gene_design_table_filename, sep='\t', index = False)
    
    
    return(nonuniq_qtl_genes_design_df)
    
############################################################################################################### 
def read_essential_genes(essential_genes_filename,
                        input_annot_gff_filename,
                        input_selected_qtls_table_filename,
                        qtl_trait_for_essential_genes,
                        include_dubious,
                        number_of_genes_to_select,
                        shuffle_genes, shuffle_genes_rand_seed = 14):
    """
    essential genes table taken from SGD: http://www.yeastgenome.org/phenotype/inviable/overview
    if qtl_trait_for_essential_genes - do not filter----
    """
    
    
    essential_genes_df = pd.read_table(essential_genes_filename, sep='\t', na_values = "None", skiprows=8)
    essential_genes_df['gene_id'] = essential_genes_df['Gene Systematic Name'].str.strip()
    
    print("%d input essential genes" % (essential_genes_df.shape[0]))
    
    # filtering relevant essential genes
    essential_genes_df = essential_genes_df[essential_genes_df['Mutant Information'] == 'null']
    essential_genes_df = essential_genes_df[essential_genes_df['Strain Background'] == 'S288C']
    
    
    print("%d essential genes after filtering SC relevant" % (essential_genes_df.shape[0]))
    
    
    # filtering genes that are not in the gff (or doubious)
    genes_gff_df = sgd_gff2dataframe(input_annot_gff_filename, ['CDS'],include_dubious=include_dubious)
    essential_genes_df = essential_genes_df[essential_genes_df['gene_id'].isin(genes_gff_df['Name'])]
    
    print("%d essential genes after filtering for not dubious" % (essential_genes_df.shape[0]))
    
    
    # filter genes out of qtls of the specified traits
    if qtl_trait_for_essential_genes is not None:
        select_qtls_df = pd.read_table(input_selected_qtls_table_filename, sep='\t', na_values = "None")
        genes_in_qtl_list = list(set(select_qtls_df['Genes Under Confidence Interval (Standard Name)'][select_qtls_df['Trait'].isin(qtl_trait_for_essential_genes)].str.cat(sep='|').split('|')))
        essential_genes_df = essential_genes_df[essential_genes_df['gene_id'].isin(genes_in_qtl_list)]
        
    print("%d essential genes selected after filter" % (essential_genes_df.shape[0]))
    
    essential_genes_df = essential_genes_df[['Gene','Gene Systematic Name', 'gene_id']]
    
    essential_genes_df.drop_duplicates(inplace=True)
    
    print("%d essential genes after removing duplicates" % (essential_genes_df.shape[0]))
    
    essential_genes_df.iloc[range(min(essential_genes_df.shape[0], number_of_genes_to_select)),:]
    
    if shuffle_genes:
        # random_state for reproducibility
        essential_genes_df = essential_genes_df.sample(frac=1,random_state=shuffle_genes_rand_seed).reset_index(drop=True)
    
    return(essential_genes_df)

    
############################################################################################################### 
def add_design_columns_to_essential_genes_old(essential_genes_df) :
    ##### join with design columns

    lookalike_essential_genes_design_df = essential_genes_df.copy()
    
    lookalike_essential_genes_design_df['set_name'] = 'essential_genes'
    lookalike_essential_genes_design_df['gene_id'] = lookalike_essential_genes_design_df['gene_id']
    lookalike_essential_genes_design_df['guide_num'] = 10
    lookalike_essential_genes_design_df['donor_mut_type'] = 'nonsense'
    lookalike_essential_genes_design_df['num_donor_variants'] = 2

    lookalike_essential_genes_design_df['mut_pos_in_guide'] = None
    lookalike_essential_genes_design_df['donor_length'] = 100
    lookalike_essential_genes_design_df['donor_seq_offsets'] = None
    lookalike_essential_genes_design_df['donor_seq_offsets'] = "[0]"
    lookalike_essential_genes_design_df['donor_seq_offsets'] = lookalike_essential_genes_design_df['donor_seq_offsets'].apply(ast.literal_eval)
    
    #for iloc in  range(lookalike_essential_genes_design_df.shape[0]):
    #    lookalike_essential_genes_design_df.set_value(iloc,'donor_seq_offsets',[-5,0,5])

    lookalike_essential_genes_design_df.to_csv(output_essential_gene_lookalike_design_table_filename, sep='\t', index = False)
    
    return(lookalike_essential_genes_design_df)



###############################################################################################################
def old_parse_essential_genes_to_lookalike_design_table(essential_genes_filename, input_annot_gff_filename,
                                                    input_selected_qtls_table_filename, number_of_genes_to_select, include_dubious, 
                                                    output_essential_gene_lookalike_design_table_filename):
    """
    # essential genes table taken from SGD: http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt
    """
    
    essential_genes_df = pd.read_table(essential_genes_filename, sep='\t', na_values = "None")
    essential_genes_df['gene_id'] = essential_genes_df['ORF_name']

    select_qtls_df = pd.read_table(input_selected_qtls_table_filename, sep='\t', na_values = "None")
    
    genes_gff_df = sgd_gff2dataframe(input_annot_gff_filename, ['CDS'],include_dubious=include_dubious)
    
    # filtering genes that are not in the gff
    essential_genes_df = essential_genes_df[essential_genes_df['gene_id'].isin(genes_gff_df['Name'])]
    

    lookalike_essential_genes_design_df = essential_genes_df
    
    lookalike_essential_genes_design_df['set_name'] = 'essential_genes'
    lookalike_essential_genes_design_df['gene_id'] = lookalike_essential_genes_design_df['gene_id']
    lookalike_essential_genes_design_df['guide_num'] = 10
    lookalike_essential_genes_design_df['donor_mut_type'] = 'nonsense'
    lookalike_essential_genes_design_df['num_donor_variants'] = 2

    lookalike_essential_genes_design_df['mut_pos_in_guide'] = None
    lookalike_essential_genes_design_df['donor_length'] = 100
    lookalike_essential_genes_design_df['donor_seq_offsets'] = None
    lookalike_essential_genes_design_df['donor_seq_offsets'] = "[0]"
    lookalike_essential_genes_design_df['donor_seq_offsets'] = lookalike_essential_genes_design_df['donor_seq_offsets'].apply(ast.literal_eval)
    
    #for iloc in  range(lookalike_essential_genes_design_df.shape[0]):
    #    lookalike_essential_genes_design_df.set_value(iloc,'donor_seq_offsets',[-5,0,5])

    lookalike_essential_genes_design_df.to_csv(output_essential_gene_lookalike_design_table_filename, sep='\t', index = False)
    
    return(lookalike_essential_genes_design_df)

###############################################################################################################
# functions for generating the oligos 
###############################################################################################################

###############################################################################################################
def generate_oligo_seq(guide_seq,donor_seq, design_df):
    """
    composing the oligos from guide and donors using the design table
    """
    ret_oligo_seq = list('N' * (1  + design_df['pos_end'].max() - design_df['pos_start'].min()))
    
    for idx,row in design_df.iterrows():
        
        if ( str(row['Name']).lower() == 'guide' ):
            ret_oligo_seq[ int(row['pos_start']) : int(row['pos_end'])+1 ] = guide_seq
        elif ( str(row['Name']).lower() == 'donor' ):
            ret_oligo_seq[ int(row['pos_start']) : int(row['pos_end'])+1 ] = donor_seq
        else:
            ret_oligo_seq[ int(row['pos_start']) : int(row['pos_end'])+1 ] = str(row['DNA_seq'])
        
    return("".join(ret_oligo_seq))

###############################################################################################################
def generate_oligo_from_guide_and_donors(input_oligo_design_table_filename,
                                         input_guide_table_filename, input_donor_table_filename,
                                         input_guide_iloc = None, input_donor_iloc = None,
                                         input_SNP_table_filename = None,
                                         output_oligos_table_filename = ''):
    
    oligo_design_df = pd.read_table(input_oligo_design_table_filename, sep='\t', na_values = "None")
    # parse oligo df

    oligo_design_df['DNA_seq'] = oligo_design_df['DNA_seq'].str.upper()

    guides_df = pd.read_table(input_guide_table_filename, sep='\t', na_values = "None")
    donors_df = pd.read_table(input_donor_table_filename, sep='\t', na_values = "None")

    print("Before filtering there are %d guides and %d donors" % (guides_df.shape[0], donors_df.shape[0]))

    # filtering the guides and donors accoring to the iloc vectors
    if input_guide_iloc is not None:
        guides_df = guides_df.iloc[input_guide_iloc]

    if input_guide_iloc is not None:
        donors_df = donors_df.iloc[input_donor_iloc]

    print("After filtering there are %d guides and %d donors" % (guides_df.shape[0], donors_df.shape[0]))

    # mergeing guides and donors on guide_id
    guide_donor_shared_columns =  [x for x in list(guides_df.columns) if x in list(donors_df.columns)]
    print('shared columns')
    print(guide_donor_shared_columns)
    oligo_df = pd.merge(guides_df, donors_df, on=guide_donor_shared_columns)
    
    # mergeing guides and donors on guide_id
    #if 'var_id' in donors_df.columns:
    #    oligo_df = pd.merge(guides_df, donors_df, on=['guide_id','var_id'])
    #elif 'Gene' in donors_df.columns:
    #    oligo_df = pd.merge(guides_df, donors_df, on=['guide_id', 'Gene'])
    #else:
    #    oligo_df = pd.merge(guides_df, donors_df, on='guide_id')

    print("joining the guides and the donors by shared columns (guide_id) creates %d oligos" % (oligo_df.shape[0]))


    # if SNP table provided add SNP details Bloom QTLs 
    if input_SNP_table_filename is not None:
        
        print("joining oligo df with SNP table: " + input_SNP_table_filename)    
        snps_df = pd.read_table(input_SNP_table_filename, sep='\t', na_values = "None")
        snps_df.rename(columns={'set_name' : 'SNP_set_name'}, inplace=True)
        #snps_df['SNP_set_name'] = snps_df['set_name']
        #snps_df.drop('set_name', axis=1, inplace=True)
        
        oligo_df = pd.merge(oligo_df, snps_df[ ['var_id', 'is_in_selected_QTL', 'Traits', 'NumTraits','SNP_set_name']  ], on = 'var_id')
        oligo_df['set_name'] = oligo_df['set_name'] + '_' + oligo_df['SNP_set_name'] 
    


    oligo_df['oligo_seq'] = ''
    oligo_df['oligo_id'] = ''
    

    for idx,row in oligo_df.iterrows():
        if (idx % 1000 == 0):
            print('parsing oligo %d out of %d' % (idx, oligo_df.shape[0]))

        oligo_df.set_value(idx,'oligo_seq', 
                           generate_oligo_seq(str(row['guide_noPAM']),
                                              str(row['donor_seq']),
                                              oligo_design_df) )

        oligo_df.set_value(idx,'oligo_id', 
                           row['set_name'] + "#" + \
                           row['guide_id'] + "#" + \
                           row['donor_id'] + "#" + \
                           "oligo_" + str(idx) )
    
    
    
    

    if output_oligos_table_filename:
        print("Saving : " + output_oligos_table_filename)
        oligo_df.to_csv(output_oligos_table_filename, sep='\t', index = False)

    return(oligo_df)


###############################################################################################################
def generate_oligo_seq_with_barcode(guide_seq, donor_seq, barcode_seq, design_df):
    """
    creates oligos from guide, donor and barcode sequences, according to design table
    """
    ret_oligo_seq = list('N' * (1  + design_df['pos_end'].max() - design_df['pos_start'].min()))
    
    for idx,row in design_df.iterrows():
        
        if ( str(row['Name']).lower() == 'guide' ):
            ret_oligo_seq[ int(row['pos_start']) : int(row['pos_end'])+1 ] = guide_seq
        elif ( str(row['Name']).lower() == 'donor' ):
            ret_oligo_seq[ int(row['pos_start']) : int(row['pos_end'])+1 ] = donor_seq
        elif ( str(row['Name']).lower() == 'barcode' ):
            ret_oligo_seq[ int(row['pos_start']) : int(row['pos_end'])+1 ] = barcode_seq
        else:
            ret_oligo_seq[ int(row['pos_start']) : int(row['pos_end'])+1 ] = str(row['DNA_seq'])
        
    return("".join(ret_oligo_seq))

###############################################################################################################
def generate_oligo_from_guide_donor_barcode(input_oligo_design_table_filename,
                                         input_guide_table_filename, input_donor_table_filename, input_barcode_table_filename = None,
                                         input_guide_iloc = None, input_donor_iloc = None,
                                         group_size = None, 
                                         output_oligos_table_filename = ''):
    """
    modified version that accepts barcode table for writing oligo
    removed input_SNP_table_filename input to streamline oligo generation
    TODO: add input_vars_pool_assignment_filename as option for pre-assigned pool-barcodes. With this file, pools are specifically assigned beforehand, without it oligos will be randomly assigned by set name
    """
    
    # read oligo design table
    oligo_design_df = pd.read_table(input_oligo_design_table_filename, sep='\t', na_values = "None")
    oligo_design_df['DNA_seq'] = oligo_design_df['DNA_seq'].str.upper()

    # read guides and donors table
    guides_df = pd.read_table(input_guide_table_filename, sep='\t', na_values = "None")
    donors_df = pd.read_table(input_donor_table_filename, sep='\t', na_values = "None")
    
    # filter guides and donors accoring to iloc vectors (if specified)
    print("Before filtering there are %d guides and %d donors" % (guides_df.shape[0], donors_df.shape[0]))
    if input_guide_iloc is not None:
        guides_df = guides_df.iloc[input_guide_iloc]
    if input_donor_iloc is not None:
        donors_df = donors_df.iloc[input_donor_iloc]
    print("After filtering there are %d guides and %d donors" % (guides_df.shape[0], donors_df.shape[0]))

    # merge guide and donor tables by guide_id
    guide_donor_shared_columns =  [x for x in list(guides_df.columns) if x in list(donors_df.columns)]
    print('shared columns')
    print(guide_donor_shared_columns)
    oligo_df = pd.merge(guides_df, donors_df, on=guide_donor_shared_columns)
    
    print("joining the guides and the donors by shared columns (guide_id) creates %d oligos" % (oligo_df.shape[0]))
    
    
    # group oligos randomly by sets specified in donor_design_table
    if input_barcode_table_filename is not None:
        # read barcodes table
        barcodes_df = pd.read_table(input_barcode_table_filename, sep=',', na_values = "None")
        # skip barcodes assigned to well 0 in table from assignment (reserved for technical use)
        barcodes_df = barcodes_df.loc[barcodes_df['Well']!=0].reset_index(drop=True)
        
        # prepare number of barcodes per group
        max_barcodes_per_well = barcodes_df.groupby('Well').size().min()
        if group_size is None:
            # use max number of barcodes for each group
            group_size = max_barcodes_per_well # smallest group of barcodes
        
        # subset barcodes table such that number of barcodes in each group matches the group size
        barcodes_df = barcodes_df.groupby('Well').head(group_size)
        barcodes_df = barcodes_df.reset_index(drop=True)

        # assign barcodes for each variant set
        # NOTE: very basic barcode assignment with randomization of oligos across each set
        # for more complex assignments, use separate script to pre-assign pool number to each oligo/variant and then fill in barcode
        to_join = []
        group_number = 0
        for group in oligo_df['set_name'].unique():
            oligo_group_df = oligo_df.loc[oligo_df['set_name']==group].sample(frac=1, random_state=1).reset_index(drop=True)
            # select available groups
            barcodes_to_assign = barcodes_df.loc[barcodes_df['Well']>group_number]
            # assign barcodes (TODO: modify to distribute oligos evenly across a defined number of wells)
            print("Assigning %d barcodes to set %s" % (len(oligo_group_df), group))
            oligo_group_df = oligo_group_df.join(barcodes_to_assign.iloc[:len(oligo_group_df)])
            # set last group number used
            group_number = oligo_group_df['Well'].max()
            to_join.append(oligo_group_df)
        
        # put guide, donor and barcode information from different sets together
        oligo_df = pd.concat(to_join).reset_index(drop=True)
        # rename barcode table columns
        oligo_df = oligo_df.rename(columns = {'Well' : 'pool', 
                                              'Barcode' : 'barcode_seq',
                                              'Unique_ID' : 'barcode_id'})
         # sort oligo table to donor table order
        oligo_df = oligo_df.set_index('donor_id').loc[donors_df['donor_id']].reset_index()[oligo_df.columns]
        
        

    else:
        print("No barcodes provided. Skipping barcode assignment.")
        oligo_df['pool'] = -1
        oligo_df['barcode_seq'] = ''
        oligo_df['barcode_id'] = 'None'
        

    # put oligo sequence together
    oligo_df['oligo_seq'] = ''
    oligo_df['oligo_id'] = ''
    
    for idx,row in oligo_df.iterrows():
        if (idx % 1000 == 0):
            print('parsing oligo %d out of %d' % (idx, oligo_df.shape[0]))

        oligo_df.set_value(idx,'oligo_seq', 
                           generate_oligo_seq_with_barcode(str(row['guide_noPAM']),
                                                           str(row['donor_seq']),
                                                           str(row['barcode_seq']),
                                                           oligo_design_df) )

        oligo_df.set_value(idx,'oligo_id', 
                           row['set_name'] + "#" + \
                           row['guide_id'] + "#" + \
                           row['donor_id'] + "#" + \
                           row['barcode_id'] + "#")# + \
                           #"oligo_" + str(idx) )

    if output_oligos_table_filename:
        print("Saving : " + output_oligos_table_filename)
        oligo_df.to_csv(output_oligos_table_filename, sep='\t', index = False)

    return(oligo_df)


###############################################################################################################
def replace_barcode_seq(seq, barcode, start, end):
    """
    takes in any sequence, searches for barcode at specified coordinates and swaps in provided barcode
    """
    seq = list(seq)
    
    if (end-start) == len(barcode):
        seq[start:end] = barcode
    else:
        warnings.warn("Specified coordinates do not match barcode length! (start {}, end {}, length {})".format(start, end, len(barcode)))
        
    return ''.join(seq)

###############################################################################################################
def assign_pool_barcode_to_oligos(input_oligo_design_table_filename,
                                  input_oligo_table_filename,
                                  input_pool_assignment_filename,
                                  input_barcode_table_filename,
                                  output_oligos_with_barcodes_filename,
                                  grouping_col = 'var_id',
                                  max_pool_size = None):
    """
    Takes oligo table and assigns pool number according to provided pool assignment table
    Randomizes oligos within each pool and assigns barcode ID and sequence, adding it to oligo sequence according
    to design specified in oligo design table (note: barcodes provided should match length in design table)
    """
    # read oligo design table
    oligo_design_df = pd.read_table(input_oligo_design_table_filename, sep='\t', na_values = "None", index_col=0)
    oligo_design_df['DNA_seq'] = oligo_design_df['DNA_seq'].str.upper()
    # get barcode coordinates from design table
    barcode_start = oligo_design_df.loc['barcode', 'pos_start']
    barcode_end = oligo_design_df.loc['barcode', 'pos_end']
    
    # read oligos table
    oligo_df = pd.read_csv(input_oligo_table_filename, sep='\t')
    
    # read pool assignment table (grouping col index must be unique)
    pool_assign_df = pd.read_csv(input_pool_assignment_filename, sep='\t', index_col=grouping_col)
    
    # read barcodes table
    barcodes_df = pd.read_csv(input_barcode_table_filename, sep=',')
    # skip barcodes assigned to well 0 in table from assignment (reserved for technical use)
    barcodes_df = barcodes_df.loc[barcodes_df['Well']!=0].reset_index(drop=True)
    
    
    # assign pool number by grouping column (variant ID by default)
    oligo_df.loc[:, 'pool'] = oligo_df[grouping_col].map(pool_assign_df['pool_num'])
    
    # assign unique barcode to each oligo by pool
    for pool, subset in oligo_df.groupby('pool'):
        # define maximum number of barcodes assignable per pool
        if max_pool_size is None:
            pool_size = len(barcodes_df.loc[barcodes_df['Well']==pool])
        else:
            pool_size = max_pool_size
            
        # check that number of oligos in subset is within pool size
        if len(subset) <= pool_size:
            # assign barcode to oligos
            oligo_df.loc[subset.index, 'barcode_seq'] = barcodes_df.loc[barcodes_df['Well']==pool, 'Barcode'][:len(subset)].values
            oligo_df.loc[subset.index, 'barcode_id'] = barcodes_df.loc[barcodes_df['Well']==pool, 'Unique_ID'][:len(subset)].values
            # add barcode sequence to oligo
            oligo_df.loc[subset.index, 'oligo_seq'] = oligo_df.loc[subset.index,:].apply(lambda x: replace_barcode_seq(x['oligo_seq'], x['barcode_seq'], barcode_start, barcode_end+1), axis=1)
            
        else:
            warnings.warn("More oligos ({}) than barcodes ({}) available for pool {}. No barcodes assigned!".format(len(subset), pool_size, pool))
    
    
    # write to output file
    if output_oligos_with_barcodes_filename:
        print("Saving : " + output_oligos_with_barcodes_filename)
        oligo_df.to_csv(output_oligos_with_barcodes_filename, sep='\t', index = False)

    return(oligo_df)

###############################################################################################################
def write_output_oligos(non_uniq_oligo_df, 
                       output_oligo_for_production_nonuniq_filename = "",
                       output_oligo_for_production_uniq_filename = "",
                       output_oligo_for_production_uniq_batch_prefix_filename = ""):
    """
    Writes final file for submission to oligonucleotide synthesis
    Arranged for matrixed oligonucleotide synthesis format
    Removed set_by_priority table (is this still useful?!)
    Oligo sequence duplication check should permit same guide-donor combinations, as long as they have a unique barcode assigned 
    """
    # write full list of oligos to file
    if output_oligo_for_production_nonuniq_filename:
        print("Saving : " + output_oligo_for_production_nonuniq_filename)
        non_uniq_oligo_df.to_csv(output_oligo_for_production_nonuniq_filename, sep='\t', index = False)

    # select unique oligo sequences
    uniq_oligo_df = non_uniq_oligo_df[['pool', 'barcode_id', 'oligo_seq', 'set_name']].drop_duplicates(subset='oligo_seq').copy()
    
    print("--- non uniq set counts:")
    print(non_uniq_oligo_df['set_name'].value_counts())

    print("--- uniq set counts:")
    print(uniq_oligo_df['set_name'].value_counts())
    
    # Sort oligos by pool and barcode numbers (barcode ID is formatted for easy pool-then-barcode numerical sort)
    uniq_oligo_df = uniq_oligo_df.sort_values(['pool', 'barcode_id']).reset_index(drop=True)
    # Get pool name by adding pool number to set name
    uniq_oligo_df['Pool name'] = 'pool_'+uniq_oligo_df['pool'].astype(str)#.str.cat(uniq_oligo_df['set_name'], sep='_')
    # add info columns for submission format
    uniq_oligo_df['Oligo number'] = np.arange(len(uniq_oligo_df))
    uniq_oligo_df['Length (bp)'] = [len(x) for x in uniq_oligo_df['oligo_seq']]
    uniq_oligo_df['Oligo name'] = uniq_oligo_df['set_name'].str.cat(uniq_oligo_df['barcode_id'], sep='_') # (set)_(well)_(barcode)
    uniq_oligo_df['Comments'] = ''
    # rename and select columns for submission format
    uniq_oligo_df = uniq_oligo_df.rename(columns={'oligo_seq':'Oligo sequence'})
    uniq_oligo_df = uniq_oligo_df[['Pool name', 'Oligo number', 'Oligo sequence', 'Length (bp)', 'Oligo name', 'Comments']]
    
    # write unique oligos to file for submission
    if output_oligo_for_production_uniq_filename:
        print("Saving : " + output_oligo_for_production_uniq_filename)
        uniq_oligo_df.to_csv(output_oligo_for_production_uniq_filename, sep='\t', index = False)
    
    return(uniq_oligo_df)


###############################################################################################################
def write_output_oligos_old2018(non_uniq_oligo_df, 
                       input_sets_by_priority_filename,
                       batch_size = 2000,
                       output_oligo_for_production_nonuniq_filename = "",
                       output_oligo_for_production_uniq_filename = "",
                       output_oligo_for_production_uniq_batch_prefix_filename = ""):
    """
    Writing the output oligos
    """
    # loading sets by priority
    sets_by_priority_df =  pd.read_table(input_sets_by_priority_filename, sep='\t', na_values = "None")
    sets_by_priority_df['set_priority'] = range(sets_by_priority_df.shape[0])
    sets_by_priority_df['oligo_num'] = 0
    
    # saving the non unique oligos
    if output_oligo_for_production_nonuniq_filename:
        print("Saving : " + output_oligo_for_production_nonuniq_filename)
        non_uniq_oligo_df.to_csv(output_oligo_for_production_nonuniq_filename, sep='\t', index = False)

    # uniq aligos
    uniq_oligo_df  = non_uniq_oligo_df[ ['oligo_seq'] ].drop_duplicates(inplace=False).copy()

    #print "# unique oligos:"
    #print uniq_oligo_df.shape[0]

    # asigning set name to the 
    uniq_oligo_df['set_name'] = "NO_SET_NAME"
    uniq_oligo_df['set_priority'] = 1000000
    
    # adding set name and priority to the oligos
    for i in range(sets_by_priority_df.shape[0]):
        set_name = str(sets_by_priority_df.ix[i, 'set_name'])
        set_priority = int(sets_by_priority_df.ix[i, 'set_priority'])
        cur_set_nonuniq_oligos = list(non_uniq_oligo_df.ix[ (non_uniq_oligo_df['set_name'] == str(set_name), 'oligo_seq')])
        uniq_oligo_df.ix[uniq_oligo_df['oligo_seq'].isin(cur_set_nonuniq_oligos), 'set_name' ] = str(set_name)
        uniq_oligo_df.ix[uniq_oligo_df['oligo_seq'].isin(cur_set_nonuniq_oligos), 'set_priority' ] = (set_priority)

    # sorting the oligos by set priorities
    uniq_oligo_df.sort_values(by=["set_priority", "set_name"], ascending = [True,True], inplace=True)
    uniq_oligo_df['oligo_production_id'] = range(uniq_oligo_df.shape[0])

    # counting each set by priority
    for i in range(sets_by_priority_df.shape[0]):
        sets_by_priority_df.ix[i, 'oligo_num' ] = sum( uniq_oligo_df['set_name'] ==  str(sets_by_priority_df.ix[i, 'set_name']) )
    
    sets_by_priority_df['oligo_num_cumsum'] =  np.cumsum(sets_by_priority_df['oligo_num'])


    print("--- non uniq set counts:")
    print(non_uniq_oligo_df['set_name'].value_counts())

    print("--- uniq set counts and cumsum (set is assigned by priority):")
    print(sets_by_priority_df)

    # saving all uniq oligos
    if output_oligo_for_production_uniq_filename:
        print("Saving : " + output_oligo_for_production_uniq_filename)
        uniq_oligo_df.to_csv(output_oligo_for_production_uniq_filename, sep='\t', index = False)

    # saving batches
    if output_oligo_for_production_uniq_batch_prefix_filename:
        batch_size = float(batch_size)
        cur_start_ind = 0
        for i in range(sets_by_priority_df.shape[0]):
            cur_cumsum = int(sets_by_priority_df.ix[i,'oligo_num_cumsum'])
            cur_set_name = str(sets_by_priority_df.ix[i,'set_name'])
            cur_end_ind = int(np.ceil(cur_cumsum / batch_size)*batch_size)
            cur_batch_filename = output_oligo_for_production_uniq_batch_prefix_filename + str(i) + "_" + cur_set_name + ".txt"
            print("Saving : " + cur_batch_filename)
            uniq_oligo_df[cur_start_ind:cur_end_ind].to_csv(cur_batch_filename, sep='\t', index = False)
            cur_start_ind = cur_end_ind
    
    return(uniq_oligo_df)
        
###############################################################################################################
def parse_oligo_id(oligo_id):
    """
    returns a dictionary of parsed oligo id
    """
    
    oligo_dict = {}
    
    substr_level1 = oligo_id.split('#')
    
    oligo_dict['set_name'] = substr_level1[0]
    oligo_dict['guide_id'] = substr_level1[1]
    oligo_dict['donor_id'] = substr_level1[2]
    oligo_dict['oligo_num_id'] = substr_level1[3]
    
    
    donor_substr = oligo_dict['donor_id'].split(':')
    
    oligo_dict['mut_type'] = donor_substr[1]
    oligo_dict['offset'] = int(donor_substr[2].replace("offset",""))
    oligo_dict['donorID'] = int(donor_substr[3].replace("donorID",""))
    oligo_dict['EditPosInGuide'] = int(donor_substr[4].replace("EditPosInGuide",""))
    
    return(oligo_dict)


###############################################################################################################
def write_small_control_oligo_set(output_oligo_for_production_nonuniq_filename, 
                                output_oligo_for_production_uniq_filename,
                                output_oligos_for_essential_genes_1bpMut_df_filename, 
                                small_control_set_guide_num,
                                small_control_set_editing_position):

    """
    selecting 1bp mutation oligo that traget a set of guides with high Azimuth score in a specific position
    1 synonymous and 1 nonsense
    """
    non_uniq_oligo_df = pd.read_table(output_oligo_for_production_nonuniq_filename, sep='\t', na_values = "None")
    uniq_oligo_df = pd.read_table(output_oligo_for_production_uniq_filename, sep='\t', na_values = "None")
    
    
    oligo_1bpMut_genes_df = pd.read_table(output_oligos_for_essential_genes_1bpMut_df_filename, sep='\t', na_values = "None")
    
    # parsing the nonuniq oligo_id
    non_uniq_1bpMut_oligo_df = non_uniq_oligo_df[ (non_uniq_oligo_df['set_name'] == 'tech_1bpMut')].copy()
    non_uniq_1bpMut_oligo_df['mut_type'] = non_uniq_1bpMut_oligo_df['oligo_id'].apply(lambda x: parse_oligo_id(x)['mut_type'])
    non_uniq_1bpMut_oligo_df['EditPosInGuide'] = non_uniq_1bpMut_oligo_df['oligo_id'].apply(lambda x: parse_oligo_id(x)['EditPosInGuide'])

    # adding azimuth score
    non_uniq_1bpMut_oligo_df = pd.merge(non_uniq_1bpMut_oligo_df, oligo_1bpMut_genes_df[['oligo_id', 'Azimuth']].drop_duplicates(), how='left', on='oligo_id')
    non_uniq_1bpMut_oligo_df[ ['guide_id','Azimuth'] ].drop_duplicates().sort_values(by=['Azimuth'], ascending=[False]). head(10)

    # selesting edit position
    non_uniq_1bpMut_oligo_df = non_uniq_1bpMut_oligo_df[non_uniq_1bpMut_oligo_df['EditPosInGuide']==small_control_set_editing_position]
    # taking specific number of guides
    non_uniq_1bpMut_oligo_df = non_uniq_1bpMut_oligo_df.ix[ non_uniq_1bpMut_oligo_df['guide_id'].isin(non_uniq_1bpMut_oligo_df[ ['guide_id','Azimuth'] ].drop_duplicates().sort_values(by=['Azimuth'], ascending=[False])[0:small_control_set_guide_num]['guide_id']),]

    # selecting one mut from each types for each guide
    syn_set = non_uniq_1bpMut_oligo_df.ix[non_uniq_1bpMut_oligo_df['mut_type'] == '1bp_synonymous_requiresOptionalNonesense','oligo_seq'][0:small_control_set_guide_num]
    nonsense_set_tmp = non_uniq_1bpMut_oligo_df.ix[non_uniq_1bpMut_oligo_df['mut_type'] == '1bp_nonsense_requiresOptionalSynonymous'  ,['guide_id','oligo_seq' ] ]
    nonsense_set = nonsense_set_tmp[ ['guide_id','oligo_seq']].groupby('guide_id').first()['oligo_seq']

    control_set_uniq_oligo_df = uniq_oligo_df[ uniq_oligo_df['oligo_seq'].isin(pd.concat([nonsense_set,syn_set])) ]


    if output_oligo_for_production_uniq_small_control_filename:
            print("Saving : " + output_oligo_for_production_uniq_small_control_filename)
            control_set_uniq_oligo_df.to_csv(output_oligo_for_production_uniq_small_control_filename, sep='\t', index = False)


###############################################################################################################
def filter_1pbMut_donor_table(input_guides_for_essential_genes_1bpMut_df_filename, input_donors_for_essential_genes_1bpMut_df_filename,
                              K_donors_for_each_guide_pos = 100, 
                              min_Azimuth_score = 0.5,
                              max_gene_pos_frac = 0.5,
                              output_donors_for_essential_genes_1bpMut_filt_df_filename = ''):
    
    """
    filtering 1bpMut donor table to have uniform number of donors for each position
    
    TODO - remove min_Azimuth_score and max_gene_pos_frac. These filters are done upstream
    """
    
    guide_df = pd.read_table(input_guides_for_essential_genes_1bpMut_df_filename, sep='\t', na_values = "None")
    donor_df = pd.read_table(input_donors_for_essential_genes_1bpMut_df_filename, sep='\t', na_values = "None")
    
    #guide_df["guide_cut_gene_pos_frac"] = guide_df["guide_cut_gene_nt_pos"] / guide_df["CDS_len_nts"]


    donor_df_synonymous = pd.merge( guide_df[ ["guide_id","guide_cut_gene_pos_frac", "Azimuth", "gene_strand"] ][ ( guide_df["guide_cut_gene_pos_frac"] <  max_gene_pos_frac).values & (guide_df["Azimuth"] > min_Azimuth_score).values ],  
                          donor_df[ (~ donor_df["donor_id"].str.contains("nonsynonymous")).values & ( donor_df["donor_id"].str.contains("_synonymous_")).values  & (~ donor_df["contain_excluded_sequences"]).values ] , on="guide_id")

    donor_df_nonsense = pd.merge( guide_df[ ["guide_id","guide_cut_gene_pos_frac", "Azimuth", "gene_strand"] ][ ( guide_df["guide_cut_gene_pos_frac"] <  max_gene_pos_frac).values & (guide_df["Azimuth"] > min_Azimuth_score).values ],  
                          donor_df[ (~ donor_df["donor_id"].str.contains("nonsynonymous")).values & ( donor_df["donor_id"].str.contains("_nonsense_")).values & (~ donor_df["contain_excluded_sequences"]).values ] , on="guide_id")

    donor_df_synonymous = donor_df_synonymous.groupby('donor_mut_pos_in_guide').head(K_donors_for_each_guide_pos)

    donor_df_synonymous_pos_uniq = donor_df_synonymous[ ['guide_id' , 'donor_mut_pos_in_guide']  ].drop_duplicates()
    donor_df_nonsense = pd.merge(donor_df_nonsense, donor_df_synonymous_pos_uniq, on = ['guide_id' , 'donor_mut_pos_in_guide'])

    donor_nonsense_and_synonymous_df = donor_df_synonymous.append(donor_df_nonsense)

    if output_donors_for_essential_genes_1bpMut_filt_df_filename:
        print("Saving : " + output_donors_for_essential_genes_1bpMut_filt_df_filename)
        donor_nonsense_and_synonymous_df.to_csv(output_donors_for_essential_genes_1bpMut_filt_df_filename, sep='\t', index = False)
    
    return(donor_nonsense_and_synonymous_df)

