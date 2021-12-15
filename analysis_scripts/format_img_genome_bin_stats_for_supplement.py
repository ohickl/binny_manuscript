import pandas as pd

bin_stats = '/Users/oskar.hickl/binny_bench/amber_data/img_genomes_all_binner_bin_stats2.tsv'
img_metadata = '/Users/oskar.hickl/binny_bench/benchamrk_assemlbies/img_data/img_genomes.tsv'

md_df = pd.read_csv(img_metadata, sep='\t', index_col=False)
md_df['Class Order'] = md_df['Class'].map(str) + ';' + md_df['Order'].map(str)
md_df['Phylum Class'] = md_df['Phylum'].map(str) + ';' + md_df['Class'].map(str)
md_df['Class Order'] = [env.split(';')[1] + ';' + env.split(';')[0] if 'Digestive system' in env
                        else env.replace(' ', ';') if 'Saline-Alkaline' in env else env for env in md_df['Class Order'].to_list()]
print(set(md_df['Class Order'].to_list()))
print(len(set(md_df['Class Order'].to_list())))

env_dict = {}

for sample, env in zip(md_df['Sample name'], md_df['Class Order']):
    if sample not in env_dict:
        env_dict[sample] = env

bin_stat_df = pd.read_csv(bin_stats, sep='\t', index_col=False)

bin_stat_df = bin_stat_df.drop(['Marker_lineage', 'N_genomes', 'N_markers', 'N_marker_sets', 'Missing_markers',
                                'Markers_1x', 'Markers_2x', 'Markers_3x', 'Markers_4x', 'Markers_>5x',
                                'Strain_heterogeneity'], axis=1)

bin_stat_df['Environment'] = [env_dict[sample] for sample in bin_stat_df['Sample']]
bin_stat_df = bin_stat_df[['Binner', 'Sample', 'Environment'] + [column for column in bin_stat_df.columns if column not in ['Binner', 'Sample', 'Environment']]]
bin_stat_df.to_csv('/Users/oskar.hickl/binny_bench/img_genomes_all_binner_bin_stats_supplement_format.tsv', sep='\t', index=False)
