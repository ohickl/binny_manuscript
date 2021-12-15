# Asses binny's read recruitment

import pysam
import glob
import pandas as pd

def load_contigs_from_fasta(fasta_path):
    contig_list = []
    with open(fasta_path, 'r') as fasta:
        for line in fasta:
            line = line.strip('\n ')
            if line.startswith('>'):
                contig_list.append(line.strip('>'))
    return contig_list


def load_fastq_read_counts(fastq_stats_path):
    with open(fastq_stats_path, 'r') as fastq_stats:
        next(fastq_stats)
        for line in fastq_stats:
            line = line.strip('\n ')
            fastq_read_count = int(line.split('\t')[3])
    return fastq_read_count


sample_list = glob.glob('binny_out/c?_t_*/[2H]*[0-9]')

read_stat_dict = {i.split('/')[1]: {} for i in sample_list}
read_stats_list = []

for sample_path in sample_list:
    sample_total_reads_mapping = 0
    sample_total_reads_mapping_single_contig = 0
    sample_data_set = sample_path.split('/')[1]
    sample_name = sample_path.split('/')[2]
    sample_bam_file_path = 'cami_data/{0}/short_read/{1}/contigs/mapping_pe_sorted.bam'.format(sample_data_set,
                                                                                               sample_name)
    sample_fastq_stats_file = 'cami_data/{0}/short_read/{1}/reads/anonymous_reads_1_stats.tsv'.format(sample_data_set,
                                                                                                      sample_name)
    sample_fastq_read_count = load_fastq_read_counts(sample_fastq_stats_file)
    sample_total_reads = sample_fastq_read_count * 2
    samfile = pysam.AlignmentFile(sample_bam_file_path, "rb")
    try:
        samfile.check_index()
    except:
        pysam.index(sample_bam_file_path)
        samfile = pysam.AlignmentFile(sample_bam_file_path, "rb")
    bin_list = glob.glob('{0}/bins_checkm/*.fasta'.format(sample_path))
    for bin in bin_list:
        bin_contigs = load_contigs_from_fasta(bin)
        for contig in bin_contigs:
            contig_read_counts = samfile.count(contig, read_callback='all')
            sample_total_reads_mapping += contig_read_counts
            if '_I' not in bin:
                sample_total_reads_mapping_single_contig += contig_read_counts
                print(sample_data_set, sample_name, bin.split('/')[-1], contig, contig_read_counts, sample_total_reads_mapping, sample_total_reads_mapping_single_contig)

    reads_mapping_ratio = round(sample_total_reads_mapping / sample_total_reads, 4)
    single_cont_reads_mapping_ratio = round(sample_total_reads_mapping_single_contig / sample_total_reads, 3)

    read_stat_dict[sample_data_set][sample_name] = [sample_total_reads_mapping, sample_total_reads_mapping_single_contig]
    read_stats_list.append([sample_data_set, sample_name, sample_total_reads_mapping, sample_total_reads_mapping_single_contig, reads_mapping_ratio, single_cont_reads_mapping_ratio, sample_total_reads])

    print(sample_data_set, sample_name, round(sample_total_reads_mapping / sample_total_reads, 3),
          round(sample_total_reads_mapping_single_contig / sample_total_reads, 3))

df_header = ['data_set', 'sample', 'reads mapping', 'single contig genome reads mapping', 'read mapping ratio', 'single contig read mapping ratio', 'sample total reads']
read_stat_df = pd.DataFrame(read_stats_list, columns=df_header)
read_stat_df = read_stat_df.sort_values(by=['data_set', 'sample'])

read_stat_df.to_csv('binny_out/binny_cami_read_rec_stats.tsv', sep='\t')

read_stat_df.describe(include='all').round({'reads mapping': 0, 'single contig genome reads mapping': 0, 'read mapping ratio': 3,
                                  'single contig read mapping ratio': 3, 'sample total reads': 0}).to_csv('binny_out/binny_cami_read_rec_stats_summary.tsv',
                                                                                                          sep='\t')
