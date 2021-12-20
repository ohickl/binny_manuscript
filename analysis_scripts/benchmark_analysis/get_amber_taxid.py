#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  28 12:22:48 2021

@author: oskar.hickl
"""

import sys


def load_sample_tax_profile(sample_tax_file):
    genome_id_to_taxid_dict = {}
    with open(sample_tax_file, 'r') as stf:
        # Skip header
        for i in range(5):
            next(stf)
        for line in stf:
            line = line.strip('\t \n').split('\t')
            if len(line) > 5:
                genome_id = line[5]
                taxid = line[0].split('.')[0]
                if genome_id in genome_id_to_taxid_dict:
                    print("Warning: {0} already in dict".format(genome_id))
                genome_id_to_taxid_dict[genome_id] = taxid
    return genome_id_to_taxid_dict


def load_contig2genome_id_dict(gsa_mapping_amber_file):
    contig2genome_id_dict = {}
    with open(gsa_mapping_amber_file, 'r') as stf:
        # Skip header
        for i in range(4):
            next(stf)
        for line in stf:
            line = line.strip('\t \n').split('\t')
            contig = line[0]
            genome_id = line[1]
            contig_length = int(line[-1])
            if contig in contig2genome_id_dict:
                print("Warning: {0} already in dict".format(contig))
            contig2genome_id_dict[contig] = [genome_id, contig_length]
    return contig2genome_id_dict


def load_contig2taxid_dict(hr_gsa_mapping_amber_file):
    contig2taxid_dict = {}
    with open(hr_gsa_mapping_amber_file, 'r') as stf:
        # Skip header
        for i in range(4):
            next(stf)
        for line in stf:
            line = line.strip('\t \n').split('\t')
            contig = line[0]
            taxid = line[2].split('.')[0]
            contig_length = int(line[-1])
            if contig in contig2taxid_dict:
                print("Warning: {0} already in dict".format(contig))
            contig2taxid_dict[contig] = [taxid, contig_length]
    return contig2taxid_dict


def write_high_res_tax_id_gs_amber_input_file(gs_amber_file, gen_id2taxid_dict, ncbi_taxid_update):
    out_file = '.'.join(gs_amber_file.split('.')[:-1]) + '_hr_tax.tsv'
    print('Writing high-res tax gs amber input to {0}.'.format(out_file))
    with open(gs_amber_file, 'r') as af:
        with open(out_file, 'w') as of:
            for i in range(4):
                header = next(af)
                of.write(header)
            for line in af:
                line = line.strip('\t \n').split('\t')
                new_taxid = gen_id2taxid_dict.get(line[1], line[2])
                line[2] = ncbi_taxid_update.get(new_taxid, new_taxid)
                of.write('\t'.join(line) + '\n')


def write_high_res_tax_id_sample_amber_input_file(amber_file, cont2largest_gen_taxid_dict, ncbi_taxid_update):   # cont2taxid_dict):
    out_file = '.'.join(amber_file.split('.')[:-1]) + '_hr_tax.tsv'
    print('Writing high-res tax amber input to {0}.'.format(out_file))
    with open(amber_file, 'r') as af:
        with open(out_file, 'w') as of:
            for i in range(4):
                header = next(af)
                if header.startswith('@@SEQUENCEID'):
                    header = header.strip('\n') + '\tTAXID\n'
                of.write(header)
            for line in af:
                line = line.strip('\t \n').split('\t')
                new_taxid = cont2largest_gen_taxid_dict.get(line[0], '2')
                line.append(ncbi_taxid_update.get(new_taxid, new_taxid))
                of.write('\t'.join(line) + '\n')


def load_bins(sample_amber_mapping):
    bin_dict = {}
    with open(sample_amber_mapping, 'r') as af:
        for i in range(4):
            header = next(af)
        for line in af:
            line = line.strip('\t \n').split('\t')
            contig = line[0]
            bin = line[1]
            if not bin in bin_dict:
                bin_dict[bin] = [contig]
            else:
                bin_dict[bin].append(contig)
    return bin_dict


def find_largest_genome_in_bins(contig2genome_id_dict, sample_tax_profile, sample_bin_dict):
    cont2largest_gen_taxid_dict = {}
    for bin, contigs in sample_bin_dict.items():
        gs_stat_dict = {}
        for contig in contigs:
            contig_genome = contig2genome_id_dict.get(contig)[0]
            contig_len = contig2genome_id_dict.get(contig)[1]
            if not contig_genome in gs_stat_dict:
                gs_stat_dict[contig_genome] = contig_len
            else:
                gs_stat_dict[contig_genome] += contig_len
        largest_genome = [None, 0]
        for k, v in gs_stat_dict.items():
            if v > largest_genome[1]:
                print('New largest genome for bin {0} is {1} with size {2}. Previous largest was {3} with size {4}.'.format(bin, k, v, largest_genome[0], largest_genome[1]))
                largest_genome = [k, v]
        try:
            largest_gen_taxid = sample_tax_profile[largest_genome[0]]
        except KeyError:
            print(largest_genome)
            print(sample_tax_profile)
            raise KeyError
        print('Largest genome for bin {0} is {1} with taxid {2}'.format(bin, largest_genome[0], largest_gen_taxid))
        for contig in contigs:
            cont2largest_gen_taxid_dict[contig] = largest_gen_taxid
    return cont2largest_gen_taxid_dict

update_taxid_dict = {'46170': '1280', '1834200': '1796646', '219334': '1423732', '1046595': '1423807', '1046597': '1423816',
                     '1122950': '1035196', '1219083': '1123383', '1240726': '1591170', '1240783': '1591174', '1243180': '1296625',
                     '1316583': '1095731', '1345697': '1921421', '335659': '1404864', '469598': '562', '654421': '265546',
                     '935543': '157910'}

write_gs = sys.argv[1]
sample_gs_mapping = sys.argv[2]
print('Gold standard mapping file is: {0}.'.format(sample_gs_mapping))
sample_tax_file = sys.argv[3]
print('Sample tax file is: {0}.'.format(sample_tax_file))
gen2taxid = load_sample_tax_profile(sample_tax_file)

if write_gs == 'True':
    print('Will write high-res tax gold standard amber file')
    write_high_res_tax_id_gs_amber_input_file(sample_gs_mapping, gen2taxid, update_taxid_dict)
else:
    sample_amber_file = sys.argv[4]
    bin_dict = load_bins(sample_amber_file)
    cont2taxid = load_contig2taxid_dict(sample_gs_mapping)
    cont2gen_id = load_contig2genome_id_dict(sample_gs_mapping)
    cont2largest_gen_taxid = find_largest_genome_in_bins(cont2gen_id, gen2taxid, bin_dict)
    write_high_res_tax_id_sample_amber_input_file(sample_amber_file, cont2largest_gen_taxid, update_taxid_dict)

