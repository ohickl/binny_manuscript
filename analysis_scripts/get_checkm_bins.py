#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 14:00:35 2021

@author: oskar.hickl
"""

import sys
from pathlib import Path
import shutil

run = Path(sys.argv[1])
sample_nr = str(run).split('_')[-1]

good_bins_list = []

checkm_out = run / 'checkm/output.txt'
sample = str(run).split('/')[-1]
if checkm_out.exists():
    with open(checkm_out, 'r') as co:
        for line in co:
            if not line.startswith(('-', '  Bin Id')):
                line = line.split()
                if int(float(line[12])) >= 70 and int(float(line[13])) <= 10:
                    good_bins_list.append(line[0])

# out_folder = Path('{0}/bins_checkm'.format(run))
out_folder = run / 'bins_checkm'
out_folder.mkdir(exist_ok=True)

second_contig = False

if len(good_bins_list) > 0:
    for bin in good_bins_list:
        bin = bin + '.fasta'
        # org_file = Path('{0}}/bins/{1}'.format(run, bin))
        org_file = run / 'bins' / bin
        dest_file = out_folder / bin
        shutil.copy(org_file, dest_file)
else:
    bin = 'dummy_bin.fasta'
    # org_file = Path('{0}}/bins/{1}'.format(run, bin))
    dest_file = out_folder / bin
    with open(dest_file, 'w') as df:
        df.write('>S{0}C0\n\n'.format(sample_nr))

