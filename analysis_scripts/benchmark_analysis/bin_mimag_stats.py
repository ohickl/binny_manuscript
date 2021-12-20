# Get bin rRNA and tRNA stats in addition ot purity and completeness to assess MIMAG quality status

import pandas as pd
import glob

# Write binner stats file
########################################################################################################################

class MAG:
    def __init__(self, name: str = '', fasta_path: str = None, gff_path: str = None) -> None:
        self.name = name
        self.contigs: dict[str, str] = {}
        if fasta_path:
            self.load_contigs_from_fasta(fasta_path)
        if gff_path:
            self.load_contig_annotations_from_gff(gff_path)
        self.stats = {'CDS': [], 'tRNA': [], 'rRNA': [], 'tmRNA': []}
        self.create_annotation_type_lists()

    def load_contigs_from_fasta(self, fasta_path) -> None:
        with open(fasta_path, 'r') as fasta:
            current_header = None
            for line in fasta:
                line = line.strip('\n ')
                # Get header
                if line.startswith('>'):
                    # Transform sub-lists of sequences from previous header to one string.
                    if current_header:
                        self.contigs[current_header]['Sequence'] = ''.join(self.contigs[current_header]['Sequence'])
                    # Set new header
                    current_header = line.strip('>')  # [3:] for vamb
                    self.contigs[current_header] = {'Sequence': []}
                # Get sequences
                else:
                    self.contigs[current_header]['Sequence'].append(line)
            # Transform sub-lists of sequences from last header to one string.
            self.contigs[current_header]['Sequence'] = ''.join(self.contigs[current_header]['Sequence'])

    def add_contig(self, identifier, sequence, overwrite=False) -> None:
        if identifier not in self.contigs:
            self.contigs[identifier] = {'Sequence': sequence}
        elif overwrite:
            print('Contig {0} already in MAG {1}. Overwriting previous entry.'.format(identifier, self.name))
            self.contigs[identifier] = {'Sequence': sequence}
        else:
            print('Contig {0} already in MAG {1}. Skipping.')

    def load_contig_annotations_from_gff(self, gff_path) -> None:
        if not self.contigs:
            print('No contigs loaded. To add contigs from a fasta file use: load_contigs_from_fasta().')
            return
        with open(gff_path, 'r') as gff:
            for line in gff:
                line = line.strip('\n ')
                line_data = line.split('\t')
                annotation_contig = line_data[0]
                if annotation_contig in self.contigs:
                    annotation_type, annotation_position, annotation_data = line_data[2], (int(line_data[3]), int(line_data[4])), line_data[-1].split(';')
                    if annotation_type not in self.contigs[annotation_contig]:
                        self.contigs[annotation_contig][annotation_type] = {}
                    if annotation_type in ['CDS', 'tRNA', 'rRNA', 'tmRNA']:
                        annotation_id = annotation_data[0].split('=')[1]
                        self.contigs[annotation_contig][annotation_type][annotation_id] = {annotation.split('=')[0]:
                                                                                               annotation.split('=')[1]
                                                                                           for annotation in annotation_data[1:]}
                    elif annotation_type == 'repeat_region':
                        crisprs_region_nr = format(len(self.contigs[annotation_contig][annotation_type].keys()), "06")
                        annotation_id = '{0}_repeat_region_{1}'.format(annotation_contig, crisprs_region_nr)
                        self.contigs[annotation_contig][annotation_type][annotation_id] = {annotation.split('=')[0]:
                                                                                               annotation.split('=')[1]
                                                                                           for annotation in annotation_data}
                    self.contigs[annotation_contig][annotation_type][annotation_id]['Position'] = annotation_position

    def create_annotation_type_lists(self) -> None:
        for contig in self.contigs:
            for annotation_type in ['CDS', 'tRNA', 'rRNA', 'tmRNA']:
                if annotation_type in self.contigs[contig]:
                    for annotation_id in self.contigs[contig][annotation_type]:
                        annotation = self.contigs[contig][annotation_type][annotation_id]['product']
                        if annotation_type == 'rRNA':
                            annotation = annotation.replace(' (partial)', '')
                        elif annotation_type == 'tRNA':
                            annotation = annotation.split('(')[0]
                        self.stats[annotation_type].append(annotation)

    def tRNAs(self) -> int:
        return len(self.stats['tRNA'])

    def unique_tRNAs(self) -> int:
        return len(set(self.stats['tRNA']))

    def tRNA_completeness(self) -> float:
        return round(len(set(self.stats['tRNA'])) / 21, 2)

    def rRNAs(self) -> int:
        return len(self.stats['rRNA'])

    def unique_rRNAs(self) -> int:
        return len(set(self.stats['rRNA']))

    def rRNA_completeness(self) -> float:
        return round(len(set(self.stats['rRNA'])) / 3, 2)

    def has_23S_rRNA(self) -> bool:
        if '23S ribosomal RNA' in self.stats['rRNA']:
            return True
        else:
            return False

    def has_5S_rRNA(self) -> bool:
        if '5S ribosomal RNA' in self.stats['rRNA']:
            return True
        else:
            return False

    def has_16S_rRNA(self) -> bool:
        if '16S ribosomal RNA' in self.stats['rRNA']:
            return True
        else:
            return False

    def has_all_rRNAs(self) -> bool:
        if '23S ribosomal RNA' in self.stats['rRNA'] and '5S ribosomal RNA' in self.stats['rRNA'] \
                and '16S ribosomal RNA' in self.stats['rRNA']:
            return True
        else:
            return False

binners = ['binny', 'concoct', 'maxbin', 'metabat', 'vamb']

for binner in binners:
    # binner_mags = glob.glob('{0}_out/img_genomes/*/bins_checkm/*.fasta'.format(binner))
    binner_mags = glob.glob('{0}_out/c?_t_*/*/bins_checkm/*.fasta'.format(binner))
    df_line_list = []
    df_header = ['data_set', 'sample', 'MAG', 'tRNAs', 'unique_tRNAs', 'tRNA_completeness', 'rRNAs', 'unique_rRNAs',
                 'rRNA_completeness', 'has_5S_rRNA', 'has_16S_rRNA', 'has_23S_rRNA']

    for mag_fasta_path in binner_mags:
        print('Working on {0}.'.format(mag_fasta_path))
        mag_data_set = mag_fasta_path.split('/')[1]
        mag_sample = mag_fasta_path.split('/')[2]
        mag_name = '.'.join(mag_fasta_path.split('/')[-1].split('.')[:-1])
        mag_gff = '/'.join(mag_fasta_path.split('/')[:3] + ['intermediary/annotation_CDS_RNA_hmms_checkm.gff']).replace(binner, 'binny')
        mag = MAG(name=mag_name, fasta_path=mag_fasta_path, gff_path=mag_gff)
        mag_line = [mag_data_set, mag_sample, mag_name, mag.tRNAs(), mag.unique_tRNAs(), mag.tRNA_completeness(),
                    mag.rRNAs(), mag.unique_rRNAs(), mag.rRNA_completeness(), mag.has_5S_rRNA(), mag.has_16S_rRNA(),
                    mag.has_23S_rRNA()]
        df_line_list.append(mag_line)

    mag_essential_RNA_stats = pd.DataFrame(df_line_list, columns=df_header)
    mag_essential_RNA_stats.to_csv('{0}_out/{0}_cami_pool_original_checkm_bins_mag_essential_RNA_stats.tsv'.format(binner), sep='\t', index=False)
    mag_essential_RNA_stats.describe(include='all').to_csv('{0}_out/{0}_cami_pool_original_checkm_bins_mag_essential_RNA_stats_summary.tsv'.format(binner),
                                                           sep='\t')

# Write combined table with AMBER and MIMAG stats for CAMI bench
########################################################################################################################

import pandas as pd

binners = ['binny', 'concoct', 'maxbin', 'metabat', 'vamb']
desc_suffix = ''
cami_out = 'cami_all'
# binners = ['binny_no_mask']
# desc_suffix = '_pool_no_mask'
# cami_out = 'cami_all_pool_no_mask'
binners_amber_names = {'binny': 'Binny', 'concoct': 'CONCOCT', 'maxbin': 'MaxBin2', 'metabat': 'MetaBAT2', 'vamb': 'VAMB'}
cami_ds_names = {'cami1_toy_high': 'c1_t_high', 'cami2_toy_airways': 'c2_t_airways', 'cami2_toy_gi': 'c2_t_gi',
                 'cami2_toy_oral': 'c2_t_oral', 'cami2_toy_skin': 'c2_t_skin', 'cami2_toy_urogen': 'c2_t_urogen'}

for binner in binners:
    binny_mimag_stats_file = '/Users/oskar.hickl/binny_bench/{0}_{1}cami_checkm_bins_mag_essential_RNA_stats.tsv'.format(binner, desc_suffix)
    amber_stats = '/Users/oskar.hickl/binny_bench/amber_data/{0}/cami_all_bin_metrics_full_bin_lineage_fix.tsv'.format(cami_out)

    min_size = 0
    max_size = 10000000
    min_pur = 0.95
    max_pur = 1.01
    min_comp = 0.9
    max_comp = 1.01
    bin_size_metric = 'seq_counts_gs'
    amber_df = pd.read_csv(amber_stats, sep='\t').query('`rank`.isnull() & `Tool` == "{0}[CheckM]" '
                                                        '& not `BINID`.str.contains("dummy_bin", na=False)'.format(binner),
                                                                                              engine='python')
    amber_df = amber_df.drop(['TAXID', 'Tool', 'filtered', 'index', 'name', 'rank', 'sample_id'], axis=1)
    amber_df['BINID'] = amber_df['BINID'].str.replace('//', '/')
    amber_df[['data_set']] = [cami_ds_names[data_set] for data_set in amber_df['data_set'].to_list()]

    amber_df[['sample']] = amber_df['BINID'].str.split('/').str[-3]
    amber_df[['MAG']] = amber_df['BINID'].str.split('/').str[-1].str.replace('.fasta', '', regex=False)
    amber_df = amber_df.drop(['BINID'], axis=1)

    mimag_df = pd.read_csv(binny_mimag_stats_file, sep='\t').query('"dummy_bin" not in `MAG`')

    mimag_df = mimag_df.sort_values(by=['data_set', 'sample', 'MAG'])
    amber_df = amber_df.sort_values(by=['data_set', 'sample', 'MAG'])

    merged_df = mimag_df.merge(amber_df, how='inner', on=['data_set', 'sample', 'MAG'])
    merged_df = merged_df.query(
                    '`precision_bp` > {0} & `recall_bp` > {1} & `precision_bp` < {2} & `recall_bp` < {3} & {4} > `{5}` > {6}'
                        .format(min_pur, min_comp, max_pur, max_comp, max_size, bin_size_metric, min_size), engine='python')

    merged_df.to_csv('/Users/oskar.hickl/binny_bench/{0}_cami{1}_checkm_bins_all_stats_NC.tsv'.format(binner, desc_suffix), sep='\t', index=False)
    merged_df.describe(include='all').to_csv('/Users/oskar.hickl/binny_bench/{0}_cami{1}_checkm_bins_all_stats_NC_summary.tsv'.format(binner, desc_suffix), sep='\t')

# Write combined table with CheckM and MIMAG stats for real-wrold bench
########################################################################################################################

import pandas as pd

binners = ['binny', 'concoct', 'maxbin', 'metabat', 'vamb']

for binner in binners:
    print(binner)
    binner_mimag_stats_file = '/Users/oskar.hickl/binny_bench/{0}_img_checkm_bins_mag_essential_RNA_stats_v03.tsv'.format(binner)
    img_stats = '/Users/oskar.hickl/binny_bench/amber_data/img_genomes_all_binner_bin_stats2.tsv'

    min_size = 0
    max_size = 10000000
    min_cont = 0
    max_cont = 5  # 0.85
    min_comp = 90
    max_comp = 101
    bin_size_metric = 'num_seqs'
    img_df = pd.read_csv(img_stats, sep='\t').query('`Binner` == "{0}"'.format(binner), engine='python')
    img_df = img_df.drop(['Marker_lineage', 'N_genomes', 'N_marker_sets', 'Missing_markers', 'Markers_1x', 'Markers_2x',
                          'Markers_3x', 'Markers_4x', 'Markers_>5x', 'Strain_heterogeneity', 'N_markers'], axis=1)

    img_df.rename(columns={'Bin': 'MAG', 'Sample': 'sample'}, inplace=True)
    img_df['MAG'] = img_df['MAG'].astype(str)

    mimag_df = pd.read_csv(binner_mimag_stats_file, sep='\t')
    mimag_df['MAG'] = mimag_df['MAG'].astype(str)

    merged_df = mimag_df.merge(img_df, how='inner', on=['sample', 'MAG'])

    merged_df = merged_df.query(
                    '`Contamination` >= {0} & `Completeness` > {1} & `Contamination` < {2} & `Completeness` <= {3} & {4} > `{5}` > {6}'
                        .format(min_cont, min_comp, max_cont, max_comp, max_size, bin_size_metric, min_size), engine='python')

    merged_df.to_csv('/Users/oskar.hickl/binny_bench/{0}_img_checkm_bins_all_stats_NC.tsv'.format(binner), sep='\t', index=False)
    merged_df.describe(include='all').to_csv('/Users/oskar.hickl/binny_bench/{0}_img_checkm_bins_all_stats_NC_summary.tsv'.format(binner), sep='\t')
