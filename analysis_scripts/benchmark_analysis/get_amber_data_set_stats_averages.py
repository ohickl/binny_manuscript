import pandas as pd

amber_run = 'cami_all_refiners_v2'  #cami_all_pool_v02  cami_all_pool_no_mask

amber_data_set_stats = '/Users/oskar.hickl/binny_bench/amber_data/{0}/results.tsv'.format(amber_run)

# Load data
df = pd.read_csv(amber_data_set_stats, sep='\t', header=0)
# Preprocess df
df = df.query('`rank`.isnull() & `Tool` != "Gold standard"', engine='python')
df = df.replace('Binny', 'binny').replace('Binny[CheckM]', 'binny[CheckM]')
df = df.drop(['rank', 'UniFrac (bp)', 'UniFrac (seq)', 'binning type'], axis=1)

# Per binner summary

frames = []
for tool in set(df['Tool'].to_list()):
    df_sum = df[df['Tool'] == tool].describe()
    df_sum.insert(loc=0, column='Tool', value=[tool for row in range(8)])
    df_sum = df_sum.rename(index={'25%': '25% percentile', '50%': '50% percentile', '75%': '75% percentile'})
    frames.append(df_sum)
pd.concat(frames, sort=False).drop(index=['count']).to_csv('/Users/oskar.hickl/binny_bench/amber_data/{0}/results_summary_v002.tsv'.format(amber_run), sep='\t')

# Per binner per data set summary
if 'pool' in amber_run:
    df = df.replace(regex=['gsa_'], value='').replace(regex=['_sample_pooled'], value='')
    df = df.rename(columns={"Sample": "data_set"})
    df = df[['data_set', 'Tool'] + [column for column in df.columns if column not in ['data_set', 'Tool']]]
    df.to_csv('/Users/oskar.hickl/binny_bench/amber_data/{0}/results_per_data_set_summary_v002.tsv'.format(amber_run), sep='\t')
else:
    df[['data_set', 'sample']] = df['Sample'].str.split("_sa", expand=True)
    df = df.replace(regex=['gsa_'], value='').replace(regex=['mple_'], value='sample_')
    frames = []
    for ds in set(df['data_set'].to_list()):
        for tool in set(df['Tool'].to_list()):
            df_sum = df.query('`data_set` == "{0}" & `Tool` == "{1}"'.format(ds, tool), engine='python').describe()
            df_sum.insert(loc=0, column='Tool', value=[tool for row in range(8)])
            df_sum.insert(loc=1, column='Data Set', value=[ds for row in range(8)])
            df_sum = df_sum.rename(index={'25%': '25% percentile', '50%': '50% percentile', '75%': '75% percentile'})
            frames.append(df_sum)
        pd.concat(frames, sort=False).drop(index=['count']).to_csv('/Users/oskar.hickl/binny_bench/amber_data/{0}/results_per_data_set_summary_v002.tsv'.format(amber_run), sep='\t')
