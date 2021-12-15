#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 14:00:43 2020

@author: oskar.hickl
"""

from pathlib import Path


def load_reduced_cami_contig_list_from_fasta(contig_file):
    remaining_contigs_list = []
    with open(contig_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                line = line.strip(' \n>')
                remaining_contigs_list.append(line)
    return remaining_contigs_list


def load_original_depth_file(org_depth_file):
    org_depth_dict = {}
    with open(org_depth_file, 'r') as f:
        # header = next(f)
        for line in f:
            line = line.strip(' \n')
            line_contig = line.split('\t')[0]
            line_depth = line.split('\t')[1]
            org_depth_dict[line_contig] = line_depth
    return org_depth_dict


def get_depth_reduction_val(org_depth_file_dict, red_depth_file):
    with open(red_depth_file, 'r') as f:
        for line in f:
            line = line.strip(' \n')
            line_contig = line.split('\t')[0]
            line_depth = line.split('\t')[1]
            line_org_depth = org_depth_file_dict.get(line_contig)
            if line_org_depth:
                depth_red_val = float(line_org_depth)/float(line_depth)
                break
    return depth_red_val


def write_more_complete_depth_file(red_contig_fasta, depth_red_val=False):
    out_path = '/'.join(str(red_contig_fasta).split('/')[:-1])
    print('Loading org depth file.')
    org_depth_file = Path(out_path) / 'anonymous_gsa_contig_depth.txt'
    org_depth_dict = load_original_depth_file(org_depth_file)
    print('Loading reduced contig list from fasta.')
    red_contig_list = load_reduced_cami_contig_list_from_fasta(red_contig_fasta)
    out_file = Path('{0}/anonymous_gsa_contig_depth_50perc_adapted.txt'.format(out_path))
    if depth_red_val:
        red_depth_file = Path('{0}/anonymous_gsa_contig_depth_50perc.txt'.format(out_path))
        depth_red_val = get_depth_reduction_val(org_depth_dict, red_depth_file)
        print('Reducing depth by factor: {0}.'.format(str(depth_red_val)))
    print('Writing output file: {0}.'.format(str(out_file)))
    with open(out_file, 'w') as of:
        for contig in red_contig_list:
            if org_depth_dict.get(contig):
                if depth_red_val:
                    of.write('\t'.join([contig] + [str(round(float(org_depth_dict.get(contig)) / depth_red_val, 10))]) + '\n')
                else:
                    of.write('\t'.join([contig] + [org_depth_dict.get(contig)]) + '\n')
            else:
                # of.write('\t'.join([contig] + 0 + '\n'))
                continue


proj_dir = Path('/scratch/users/ohickl/binning/binny_eval/cami_data/short_read')
fasta_files = list(proj_dir.glob('2017.12.04_18.45.54_sample_*/contigs/anonymous_gsa_50perc.fasta'))

for fasta in fasta_files:
    sample = str(fasta).split('/')[-3]
    print('Working on {0}.'.format(sample))
    write_more_complete_depth_file(fasta, depth_red_val=False)

