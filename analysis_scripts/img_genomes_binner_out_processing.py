
img_checkm_out_txt = '/Users/oskar.hickl/binny_bench/amber_data/img_genomes_all_binner_checkm_out.txt'
img_checkm_out_tsv = '/Users/oskar.hickl/binny_bench/amber_data/img_genomes_all_binner_bin_stats_all_bins.tsv'
img_bin_fasta_stats_tsv = '/Users/oskar.hickl/binny_bench/amber_data/img_genomes_all_binner_bin_fasta_stats.tsv'

bin_fastas_dict = {}
with open(img_bin_fasta_stats_tsv, 'r') as f:
    bin_fasta_stats_header = next(f)
    bin_fasta_stats_header = bin_fasta_stats_header.strip('\t ').split('\t')
    bin_fasta_stats_header = bin_fasta_stats_header[3:11] + [bin_fasta_stats_header[-3]]
    for line in f:
        line = line.strip('\t ').split('\t')
        line_binner = line[0].split('/')[0].split('_')[0]
        line_sample = line[0].split('/')[2]
        line_bin = line[0].split('/')[-1].replace('.fasta', '')
        bin_fastas_dict[';'.join([line_binner, line_sample, line_bin])] = line[3:11] + [line[-3]]

header = '\t'.join(['Binner', 'Sample', 'Bin', 'Marker_lineage', 'N_genomes', 'N_markers', 'N_marker_sets', 'Missing_markers',
                    'Markers_1x', 'Markers_2x', 'Markers_3x', 'Markers_4x', 'Markers_>5x', 'Completeness', 'Contamination',
                    'Strain_heterogeneity'] + bin_fasta_stats_header) + '\n'

with open(img_checkm_out_txt, 'r') as f:
    with open(img_checkm_out_tsv, 'w') as of:
        of.write(header)
        for line in f:
            line = line.strip('\t ').split()
            line_sample = line[0].split('/')[2]
            line_binner = line[0].split('/')[0].split('_')[0]
            line_bin = line[1]
            if ';'.join([line_binner, line_sample, line_bin]) in bin_fastas_dict:
                line = '\t'.join([line_binner] + [line_sample] + [line_bin] +
                                 ['_'.join([line[2], line[3]])] + line[4:] +
                                 bin_fastas_dict[';'.join([line_binner, line_sample, line_bin])]) + '\n'
                of.write(line)
