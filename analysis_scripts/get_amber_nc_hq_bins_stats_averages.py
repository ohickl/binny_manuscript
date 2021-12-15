import pandas as pd

amber_data_set_stats = '/Users/oskar.hickl/binny_bench/amber_data/cami_all/bin_metrics.tsv'

min_size = 0
max_size = 10000000
min_pur = 0.9
max_pur = 1.01  # 0.85
min_comp = 0.7
max_comp = 1.01
rank = 'NaN'
bin_size_metric = 'seq_counts_gs'

# Load data
df = pd.read_csv(amber_data_set_stats, sep='\t', header=0)
# Preprocess df
df = df.query('`precision_bp` > {0} & `recall_bp` > {1} & `precision_bp` < {2} & `recall_bp` < {3} & {4} > `{5}` > {6}'
                ' & `rank`.isnull() & `Tool` != "Gold standard"'.format(min_pur, min_comp, max_pur, max_comp, max_size,
                                                                        bin_size_metric, min_size), engine='python')

df = df.replace('Binny', 'binny').replace('Binny[CheckM]', 'binny[CheckM]')

df = df.drop(['rank', 'TAXID', 'filtered', 'index', 'name'], axis=1)

frames = []
for tool in set(df['Tool'].to_list()):
    for quality in ['HQ', 'NC']:
        if quality == 'NC':
            min_pur = 0.95
            min_comp = 0.9
        else:
            min_pur = 0.9
            min_comp = 0.7
        df_sum = df[df['Tool'] == tool].query('`precision_bp` > {0} & `recall_bp` > {1}'.format(min_pur, min_comp)).describe()
        df_sum.insert(loc=0, column='Tool', value=[tool for row in range(8)])
        df_sum.insert(loc=1, column='Quality', value=[quality for row in range(8)])
        frames.append(df_sum)
pd.concat(frames, sort=False).to_csv('/Users/oskar.hickl/binny_bench/amber_data/cami_all/hq_nc_bin_metrics_summary_v001.tsv', sep='\t')
