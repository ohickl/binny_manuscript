#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Fri Apr  3 08:22:12 2020

@author: oskar.hickl
'''

import sys
from pathlib import Path


def load_fasta(fasta):
    """

    Parameters
    ----------
    fasta : PosixPath/STR
        Path to input fasta.
    Returns
    -------
    Dict with headers as keys and sequences as list.

    """
    sequence_dict = {}
    with open(fasta) as f:
        current_header = None
        for line in f:
            # Get header
            if line.startswith('>'):
                # Transform sub lists of sequences from previous header to one list.
                if current_header:
                    sequence_dict[current_header] = ''.join(sequence_dict[current_header])
                line = line.strip('\n ')
                line = line.strip('>')
                sequence_dict[line] = []
                current_header = line
            # Get sequences
            else:
                line = line.strip('\n ')
                sequence_dict[current_header].append(line)
        # Transform sub lists of sequences from last header to one list.
        sequence_dict[current_header] = ''.join(sequence_dict[current_header])
    return sequence_dict


def get_cami_fastas(assembly_dict, gsa_mapping, out_dir):
    # Open gsa mapping.
    with open(gsa_mapping) as f:
        header = next(f)
        # For each line get contig and genome.
        for line in f:
            line = line.strip('\n \t')
            contig = line.split()[0]
            genome = line.split()[1]
            genome_file_name = 'gsa_' + genome + '.fasta'
            out_genome = out_dir / genome_file_name
            with open(out_genome, 'a') as out_file:
                out_file.write('>' + contig + '\n' + assembly_dict[contig] + '\n')
    print('Written gsa genome fastas.') 


input_assembly = Path(sys.argv[1])
input_gsa_mapping = Path(sys.argv[2])
output_dir = Path(sys.argv[3])

input_assembly_dict = load_fasta(input_assembly)
get_cami_fastas(input_assembly_dict, input_gsa_mapping, output_dir)

