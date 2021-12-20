#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 14:49:35 2020

@author: oskar.hickl
"""

import sys
from pathlib import Path

gs_map_in = Path(sys.argv[1])
gs_map_amber_format = Path(sys.argv[2])
sample = str(gs_map_in).split("/")[2]
sample = sample.split("_")[2]+"_"+sample.split("_")[3]

with open(gs_map_amber_format,'w+') as o:
    o.write("@Version:0.9.1\n"+"@SampleID:gsa_cami2_toy_gut_"+sample+"\n\n"+"@@SEQUENCEID\tBINID\tTAXID\t_LENGTH\n")

with open(gs_map_in) as i:
    next(i)
    for line in i:
        line=line.rstrip()
        contig = line.split()[0]
        bin = line.split()[1]
        taxid = line.split()[2]
        start = int(line.split()[5])
        end = int(line.split()[6])
        length = end-start
        with open(gs_map_amber_format,'a') as o:
            o.write(contig+"\t"+bin+"\t"+taxid+"\t"+str(length)+"\n")
