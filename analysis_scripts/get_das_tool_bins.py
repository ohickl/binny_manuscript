#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 08:22:12 2020

@author: oskar.hickl
"""

import sys
from pathlib import Path

assembly=Path(sys.argv[1])
binny_out=Path(sys.argv[2])
out_dir=Path(sys.argv[3])

assembly_dict={}

def load_fasta(FASTA,DICT):
    """
    
    Parameters
    ----------
    FASTA : STR
        Path to input fasta.
    DICT : Dict
       Dict to put fasta headers as keys and sequences as list into.

    Returns
    -------
    Dict with headers as keys and sequences as list.

    """
    with open(FASTA) as f:
        current_header=None
        for line in f:
            # Get header
            if line.startswith(">"):
                # Transform sublists of sequences from previous header to one list.
                if current_header:
                    DICT[current_header]  = "".join(DICT[current_header])
                line=line.rstrip()
                line=line.strip(">")
                DICT[line]=[]
                current_header=line
            # Get sequences
            else:
                line=line.rstrip()
                DICT[current_header].append(line)
        # Transform sublists of sequences from last header to one list.
        DICT[current_header] = "".join(DICT[current_header])
    return DICT

load_fasta(assembly,assembly_dict)

# Open binny output.
with open(binny_out) as b:
    # For each line get contig and bin it was put into.
    for line in b:
        line = line.rstrip()
        contig = line.split()[0]
        bin = line.split()[1]
        # This needs fixing. Should dynamically select bin prefix.
        bin_file_name = "das_tool_with_cc."+bin+".fasta"
        out_bin = out_dir / bin_file_name
        out = open(out_bin,"a")
        out.write('>'+contig+"\n"+assembly_dict[contig]+'\n')
