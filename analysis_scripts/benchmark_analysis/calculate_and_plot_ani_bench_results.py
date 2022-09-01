
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import glob
from collections import defaultdict


def get_genomes_with_high_closest_ani(path_to_ani_outs):
    ani_file_paths = glob.glob(f'{path_to_ani_outs}/*/short_read/*/contigs/genomes/genomes_ANI.tsv')
    genome2genomes_ani_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for path in ani_file_paths:
        path_comps = path.split('/')
        data_set = path_comps[-6]
        sample = '_'.join(path_comps[-4].split('_')[-2:]).replace('H_S00', 'sample_')
        with open(path, 'r') as f:
            for line in f:
                line = line.rstrip().split('\t')
                query_genome = line[0].split('/')[-1].replace('.fasta', '')
                ref_genome = line[1].split('/')[-1].replace('.fasta', '')
                line_ani = float(line[2])
                if query_genome != ref_genome:
                    genome2genomes_ani_dict[data_set][sample][query_genome].append((ref_genome, line_ani))
        for sample, sample_genomes in genome2genomes_ani_dict[data_set].items():
            for genome, genome_ani_list in sample_genomes.items():
                genome_ani_list.sort(key=lambda y: y[1], reverse=True)
    return genome2genomes_ani_dict


def load_binner_nc_hq_bins(cami_amber_out_bin_stats):
    line_count = 0
    binner_bin_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    with open(cami_amber_out_bin_stats, 'r') as f:
        header = next(f)
        for line in f:
            line = line.rstrip().split('\t')
            if not line[10]:
                line_count += 1
                binner = line[2]
                genome = line[4]
                purity = float(line[8])
                completeness = float(line[11])
                n_seqs_gs = int(line[13])
                n_seqs_tp = int(line[18])
                n_seqs_used = n_seqs_gs
                sample = '_'.join(line[19].split('_')[-2:])
                data_set = line[20].replace('toy', 't').replace('cami', 'c')
                if purity > 0.95 and completeness > 0.90:
                    binner_bin_stats[binner][data_set][sample]['NC'].append((genome, n_seqs_used))
                    binner_bin_stats[binner][data_set][sample]['HQ'].append((genome, n_seqs_used))
                elif purity > 0.90 and completeness > 0.70:
                    binner_bin_stats[binner][data_set][sample]['HQ'].append((genome, n_seqs_used))
    return binner_bin_stats


def get_binner_high_close_rel_ani_bins(g2g_ani_dict, binner_amber_stat_dict, genome_min_contig_sizes=(0, 5)):
    binner_high_close_rel_ani_bins_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
    for binner, data_sets in binner_amber_stat_dict.items():
        for data_set, samples in data_sets.items():
            for sample, qual_trheshs in samples.items():
                for qual, genome_list in qual_trheshs.items():
                    for genome in genome_list:
                        gen_name = genome[0]
                        if g2g_ani_dict[data_set][sample][gen_name]:
                            for ani in np.arange(90, 100.1, 0.1):
                                ani = round(ani, 1)
                                for genome_min_contig_size in genome_min_contig_sizes:
                                    if g2g_ani_dict[data_set][sample][gen_name][0][1] > ani and genome[1] >= genome_min_contig_size:
                                        binner_high_close_rel_ani_bins_dict[binner][qual][genome_min_contig_size][ani] += 1
        for qual in ['NC', 'HQ']:
            for ani in np.arange(90, 100.1, 0.1):
                ani = round(ani, 1)
                for genome_min_contig_size in genome_min_contig_sizes:
                    if not binner_high_close_rel_ani_bins_dict[binner][qual][genome_min_contig_size][ani]:
                        binner_high_close_rel_ani_bins_dict[binner][qual][genome_min_contig_size][ani] = 0
    return binner_high_close_rel_ani_bins_dict


# Path to cami_all_ani_bench_input.zip data (obtainable from zenodo)
fastani_path = '/Users/oskar.hickl/binny_bench/fastani'
cami_amber_bin_stats_processed = f'{fastani_path}/cami_all_samples_all_bins_stats.tsv'

g2g_ani_dict = get_genomes_with_high_closest_ani(fastani_path)
binner_bin_stat_dict = load_binner_nc_hq_bins(cami_amber_bin_stats_processed)

min_cont_sizes = [0, 6]

binner_high_ani_neigh_stats = get_binner_high_close_rel_ani_bins(g2g_ani_dict, binner_bin_stat_dict,
                                                                 genome_min_contig_sizes=min_cont_sizes)

sns.set_theme()

for cont_size in min_cont_sizes:
    data = [[binner, ani, mag_count, qual, mint_cont_size] for binner in binner_high_ani_neigh_stats.keys()
                                        for qual in ['HQ', 'NC']
                                        for mint_cont_size in min_cont_sizes
                                        for ani, mag_count in binner_high_ani_neigh_stats[binner][qual][mint_cont_size].items()
                                        if (']' in binner) and ani >= 90.0 and mint_cont_size == cont_size]  #  and ani >= 99.8   or binner == 'Gold standard'


    df = pd.DataFrame(data, columns=['Binning Method', 'Minimum ANI to the most similar genome in the same sample', 'Number of MAGs',
                                     'Quality', 'Minimum Genome Fragmentation [#contigs]'])


    palette = {'binny[CheckM]': '#afafe9', 'CONCOCT[CheckM]': '#afe9af', 'MaxBin2[CheckM]': '#ffeeaa', 'MetaBAT2[CheckM]': '#ff8080',
               'MetaDecoder[CheckM]': '#B06330', 'SemiBin[CheckM]': '#30B097', 'VAMB[CheckM]': '#e580ff', 'Gold standard': '#808080'}

    ax = sns.relplot(data=df, x="Minimum ANI to the most similar genome in the same sample", y="Number of MAGs", col="Quality",
                row="Minimum Genome Fragmentation [#contigs]", hue="Binning Method",
                ci=None, kind="line", facet_kws={"sharey": False}, palette=palette)

    plt.show()
    # plt.savefig(f'/path/to/out/file_{cont_size}.pdf', bbox_inches='tight',
    #                                 edgecolor='black')
