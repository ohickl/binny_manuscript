# Create binner recall violin plots for CAMI or IMG benchs

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
from matplotlib.patches import PathPatch
import seaborn as sns
import numpy as np
from upsetplot import UpSet
import re


sns.set_theme()
palette = 'deep'

img_checkm_out_tsv = '/Users/oskar.hickl/binny_bench/amber_data/img_genomes_all_binner_bin_stats2.tsv'
img_metadata = '/Users/oskar.hickl/binny_bench/benchamrk_assemlbies/img_data/genomecart_export_metabat_img_list.xls'
img_assembly_data = '/Users/oskar.hickl/binny_bench/amber_data/img_genomes_assembly_stats.tsv'
cami_all_bin_tsv = '/Users/oskar.hickl/binny_bench/amber_data/cami_all/bin_metrics.tsv'
cami_all_sample_tsv = '/Users/oskar.hickl/binny_bench/amber_data/cami_all/results.tsv'

# Choose either CAMI or img data
data_set_used = img_checkm_out_tsv  # img_checkm_out_tsv cami_all_bin_tsv

min_size = 0
max_size = 1000000

if data_set_used == img_checkm_out_tsv:
    ds_list = ['Plants:\nRoots', 'Terrestrial:\nDeep subsurface', 'Bioreactor:\nAnaerobic',
               'Digestive system:\nArthropoda',
               'Aquatic:\nFreshwater', 'Digestive system:\nAnnelida', 'Terrestrial:\nSoil',
               'Porifera:\nUnclassified', 'Lab enrichment:\nDefined media',
               'Aquatic:\nThermal springs', 'Wastewater:\nNutrient removal', 'Aquatic:\nMarine',
               'Aquatic:\nNon-marine\nSaline-Alkaline',
               'Wastewater:\nAnaerobic digestor', 'Solid waste:\nFeedstock', 'Wastewater:\nActivated Sludge',
               'Digestive system:\nMammals']
else:
    ds_list = ['cami1_toy_high', 'cami2_toy_airways', 'cami2_toy_gi', 'cami2_toy_oral', 'cami2_toy_skin',
               'cami2_toy_urogen']

# for data_set_to_plot in ds_list:  ### For single data set plotting, put everything belwo in for loop
df = pd.read_csv(data_set_used, sep='\t', index_col=False)
if data_set_used == img_checkm_out_tsv:
    max_cont = 5
    min_cont = 1.01
    min_comp = 90
    max_comp = 1.01
    binner_col_name = 'Binner'
    sample_col_name = 'Sample'
    n_seq_col_name = 'num_seqs'
    query_string = '`Contamination` < {0} & `Completeness` > {1} & `Binner` == "{2}" & `Sample` == "{3}"'
    pur_quals, comp_quals = [5, 10], [90, 70]
    md_df = pd.read_csv(img_metadata, sep='\t', index_col=False)
    md_df['Class Order'] = md_df['Class'].map(str) + ':\n' + md_df['Order'].map(str)
    md_df['Phylum Class'] = md_df['Phylum'].map(str) + ' ' + md_df['Class'].map(str)
    md_df['Class Order'] = [env.split(':\n')[1] + ':\n' + env.split(':\n')[0] if 'Digestive system' in env
                            else env.replace(' ', '\n') if 'Saline-Alkaline' in env else env for env in md_df['Class Order'].to_list()]
    ######### For individual data set plots
    # md_df = md_df.query('`Class Order`.str.contains("{0}")'
    #                     ' & `Class Order`.str.contains("{1}")'.format(data_set_to_plot.split('\n')[0],
    #                                                                   data_set_to_plot.split('\n')[1]),
    #                                                                   engine='python')
    # ds_sample_list = list(str(i) for i in set(md_df['Sample name'].to_list()) if str(i) != 'nan')
    #########

    assembly_stat_df = pd.read_csv(img_assembly_data, sep='\t', index_col=False)
    assembly_stat_df['sample'] = [file_path[1] for file_path in assembly_stat_df['file'].str.split('/').to_list()]

    ######### For individual data set plots
    # assembly_stat_df = assembly_stat_df.query('`sample`.isin({0})'.format(ds_sample_list), engine='python')
    # df = df.query('`Sample`.isin({0})'.format(ds_sample_list), engine='python')
    #########

    img_sample_data_set_size_stats = {'Assemblies': {}}
    for sample in set(assembly_stat_df['sample'].to_list()):
        img_sample_data_set_size_stats['Assemblies'][sample] = assembly_stat_df.query('`sample` == "{0}" '.format(sample),
                                                                                      engine='python')['sum_len'].to_list()[0]
else:
    max_cont = 0.95
    min_cont = 1.01
    min_comp = 0.90
    max_comp = 1.01
    binner_col_name = 'Tool'
    sample_col_name = 'sample_id'
    n_seq_col_name = 'seq_counts_gs'
    query_string = '`precision_bp` > {0} & `recall_bp` > {1} & `Tool` == "{2}" & `sample_id` == "{3}" & `rank`.isnull()'
    pur_quals, comp_quals = [0.95, 0.90], [0.90, 0.70]
    sample_stat_df = pd.read_csv(cami_all_sample_tsv, sep='\t', index_col=False)

    ######### For single data set plotting
    # sample_stat_df = sample_stat_df.query('`Sample`.str.contains("{0}")'.format(data_set_to_plot), engine='python')
    # df = df.query('`sample_id`.str.contains("{0}")'.format(data_set_to_plot), engine='python')
    #########

    cami_sample_data_set_size_stats = {'Gold standard': {}}
    for sample in set(df.sample_id.to_list()):
        cami_sample_data_set_size_stats['Gold standard'][sample] = sum(df.query('`Tool` == "Gold standard" '
                                                                                '& `sample_id` == "{0}" '
                                                                                '& `rank`.isnull() '
                                                                                '& `{1}` >= {2}'
                                                                                '& `{1}` <= {3}'
                                                                                .format(sample, n_seq_col_name,
                                                                                        min_size, max_size),
                                                                                engine='python')['length_gs'].to_list())

    df['data_set'] = ['_'.join(sample.split('_')[1:4]) for sample in df.sample_id.to_list()]
    for data_set in set(df.data_set.to_list()):
        cami_sample_data_set_size_stats['Gold standard'][data_set] = sum(df.query('`Tool` == "Gold standard" '
                                                                                  '& `data_set` == "{0}" '
                                                                                  '& `rank`.isnull()'
                                                                                  '& `{1}` >= {2}'
                                                                                  '& `{1}` <= {3}'
                                                                                  .format(data_set, n_seq_col_name,
                                                                                          min_size, max_size),
                                                    engine='python')['length_gs'].to_list())



binner_col = []
sample_col = []
bin_count_col = []
qual_col = []
data_set_col = []
perc_binned_col = []

binner_bin_count_dict = {binner: {'NC': 0, 'HQ': 0} for binner in set(df[binner_col_name].to_list())
                         if ("Gold standard" in binner or "[CheckM]" not in binner) and "Gold" not in binner}

for binner in binner_bin_count_dict.keys():
    for sample in set(df[sample_col_name].to_list()):
        for max_cont, min_comp in zip(pur_quals, comp_quals):
            if data_set_used == cami_all_bin_tsv and binner != "Gold standard":
                binner_query = binner + '[CheckM]'
            else:
                binner_query = binner
            binner_df = df.query(query_string.format(max_cont, min_comp, binner_query, sample), engine='python')
            if data_set_used == cami_all_bin_tsv:
                binner_df['Tool'] = [binner[0] for binner in binner_df['Tool'].str.split('[').to_list()]

            if min_size > 0:
                binner_df = binner_df.query('`seq_counts_gs` >= {0}'.format(min_size), engine='python')

            if max_size < 1000000:
                binner_df = binner_df.query('`seq_counts_gs` <= {0}'.format(max_size), engine='python')
            binner_col.append(binner)
            sample_col.append(sample)
            if data_set_used == img_checkm_out_tsv:
                data_set_col.append(md_df[md_df['Sample name'] == sample]['Class Order'].to_list()[0])
            else:
                data_set_col.append('_'.join(sample.split('_')[1:4]))
            bin_count_col.append(len(binner_df.index))
            if max_cont == pur_quals[0] and min_comp == comp_quals[0]:
                bin_qual = 'NC'
                qual_col.append(bin_qual)
            elif max_cont == pur_quals[1] and min_comp == comp_quals[1]:
                bin_qual = 'HQ'
                qual_col.append(bin_qual)
            binner_bin_count_dict[binner][bin_qual] += len(binner_df.index)

            if data_set_used == cami_all_bin_tsv:
                if binner not in cami_sample_data_set_size_stats:
                    cami_sample_data_set_size_stats[binner] = {'NC': {}, 'HQ': {}}
                for sample in set(binner_df.sample_id.to_list()):
                    if not binner_df.query('`Tool` == "{0}" & `sample_id` == "{1}" & `rank`.isnull()'.format(binner, sample),
                                 engine='python')['total_length'].to_list():
                        raise Exception
                    cami_sample_data_set_size_stats[binner][bin_qual][sample] = sum(
                        binner_df.query('`Tool` == "{0}" & `sample_id` == "{1}" & `rank`.isnull()'.format(binner, sample),
                                 engine='python')['total_length'].to_list())

                binner_df['data_set'] = ['_'.join(sample.split('_')[1:4]) for sample in binner_df.sample_id.to_list()]
                for data_set in set(binner_df.data_set.to_list()):
                    if data_set not in cami_sample_data_set_size_stats[binner][bin_qual]:
                        cami_sample_data_set_size_stats[binner][bin_qual][data_set] = sum(
                            binner_df.query('`Tool` == "{0}" & `data_set` == "{1}" & `rank`.isnull()'.format(binner, data_set),
                                     engine='python')['total_length'].to_list())
                    else:
                        cami_sample_data_set_size_stats[binner][bin_qual][data_set] += sum(
                            binner_df.query('`Tool` == "{0}" & `data_set` == "{1}" & `rank`.isnull()'.format(binner, data_set),
                                     engine='python')['total_length'].to_list())
            else:
                if binner not in img_sample_data_set_size_stats:
                    img_sample_data_set_size_stats[binner] = {'NC': {}, 'HQ': {}}
                img_sample_data_set_size_stats[binner][bin_qual][sample] = sum(
                    binner_df.query('`Binner` == "{0}" & `Sample` == "{1}"'.format(binner, sample),
                                    engine='python')['sum_len'].to_list())

if data_set_used == img_checkm_out_tsv:
    for binner in img_sample_data_set_size_stats:
        if binner != 'Assemblies':
            for sample in set(assembly_stat_df['sample'].to_list()):
                if sample not in img_sample_data_set_size_stats[binner]['HQ']:

                    img_sample_data_set_size_stats[binner]['HQ'][sample] = 0
                if sample not in img_sample_data_set_size_stats[binner]['NC']:
                    img_sample_data_set_size_stats[binner]['NC'][sample] = 0


if data_set_used == img_checkm_out_tsv:
    stat_dict_used = img_sample_data_set_size_stats
    reference_name = 'Assemblies'
else:
    stat_dict_used = cami_sample_data_set_size_stats
    reference_name = 'Gold standard'

for i in stat_dict_used:
    if i != 'Assemblies' and i != 'Gold standard':
        print(i, len(stat_dict_used[i]['NC']))

for binner in binner_bin_count_dict:
    for bin_qual in ['NC', 'HQ']:
        percs_recovered_sample = [stat_dict_used[binner][bin_qual].get(sample, 0) /
                                  (stat_dict_used[reference_name][sample] + 1e-10) for sample in
                                  stat_dict_used[reference_name].keys() if sample.startswith('gsa_') or data_set_used == img_checkm_out_tsv]
        percs_recovered_sample_average = round(np.average(percs_recovered_sample) * 100, 1)
        percs_recovered_sample_sem = round(np.std(percs_recovered_sample, ddof=1) / np.sqrt(np.size(percs_recovered_sample)) * 100, 1)
        stat_dict_used[binner][bin_qual]['Average Recall per Sample [bp]'] = [percs_recovered_sample_average,
                                                                                               percs_recovered_sample_sem]
        if data_set_used == cami_all_bin_tsv:
            percs_recovered_data_set = [stat_dict_used[binner][bin_qual].get(sample, 0) /
                                        (stat_dict_used[reference_name][sample] + 1e-10) for sample in
                                      stat_dict_used[reference_name].keys() if not sample.startswith('gsa_')]
            percs_recovered_data_set_average = round(np.average(percs_recovered_data_set) * 100, 1)
            percs_recovered_data_set_sem = round(np.std(percs_recovered_data_set, ddof=1) / np.sqrt(np.size(percs_recovered_data_set)) * 100, 1)

            stat_dict_used[binner][bin_qual]['Average Recall per Data Set [bp]'] = [percs_recovered_data_set_average,
                                                                                                     percs_recovered_data_set_sem]

for k in stat_dict_used.keys():
    if k != 'Gold standard' and k != 'Assemblies':
        if data_set_used != img_checkm_out_tsv:
            print(k, 'Data Set Recall', stat_dict_used[k]['NC']['Average Recall per Data Set [bp]'])
        print(k, 'Sample Recall NC:', stat_dict_used[k]['NC']['Average Recall per Sample [bp]'],
              'Sample Recall HQ:', stat_dict_used[k]['HQ']['Average Recall per Sample [bp]'])

plot_data_dict = {'Binner': binner_col, 'Sample': sample_col, 'Bin Count': bin_count_col, 'Bin Quality': qual_col,
                  'Data Set': data_set_col}

if data_set_used == img_checkm_out_tsv:
    missing_binner_list = ['concoct', 'binny', 'maxbin', 'metabat', 'vamb']
else:
    missing_binner_list = ['CONCOCT', 'Binny', 'MaxBin2', 'MetaBAT2', 'VAMB']

for missing_binner in missing_binner_list:
    if missing_binner not in binner_col:
        tmp_list = set([ds + '|' + dd for ds, dd in zip(sample_col, data_set_col)])
        tmp_sample_col = [i.split('|')[0] for i in tmp_list]
        tmp_ds_col = [i.split('|')[1] for i in tmp_list]
        for dict_sample, dict_ds in zip(tmp_sample_col, tmp_ds_col):
            for bin_qual in ['NC', 'HQ']:
                plot_data_dict['Binner'].append(missing_binner)
                plot_data_dict['Sample'].append(dict_sample)
                plot_data_dict['Bin Count'].append(0)
                plot_data_dict['Bin Quality'].append(bin_qual)
                plot_data_dict['Data Set'].append(dict_ds)
                stat_dict_used[missing_binner] = {'NC': {dict_sample: 0,
                                                        'Average Recall per Sample [bp]': [0, np.nan]},
                                                 'HQ': {dict_sample: 0,
                                                        'Average Recall per Sample [bp]': [0, np.nan]}}


plot_data_dict['Recall'] = [stat_dict_used[binner][bin_qual].get(sample, 0) /
                            (stat_dict_used[reference_name][sample] + 1e-10) for binner, sample, bin_qual in
                                  zip(plot_data_dict['Binner'], plot_data_dict['Sample'], plot_data_dict['Bin Quality'])
                            if sample.startswith('gsa_') or data_set_used == img_checkm_out_tsv]

df2 = pd.DataFrame(plot_data_dict)
df2 = df2.sort_values(by=['Binner'])

# Comment for stat dfs
df2['Binner'] = df2['Binner'].map(str) + '_' + df2['Bin Quality'].map(str)
# To create stats dfs
########################################################################################################################
# df_stats = df2.replace('Binny', 'binny').replace('Binny[CheckM]', 'binny[CheckM]').replace('concoct', 'CONCOCT')\
#               .replace('maxbin', 'MaxBin2').replace('metabat', 'MetaBAT2').replace('vamb', 'VAMB')
# frames = []
# for tool in set(df_stats['Binner'].to_list()):
#     for quality in ['HQ', 'NC']:
#         df_sum = df_stats[df_stats['Binner'] == tool].query('`Bin Quality` == "{0}"'.format(quality)).describe()
#         df_sum.insert(loc=0, column='Tool', value=[tool for row in range(8)])
#         df_sum.insert(loc=1, column='Quality', value=[quality for row in range(8)])
#         df_sum = df_sum.rename(index={'25%': '25% percentile', '50%': '50% percentile', '75%': '75% percentile'})
#         frames.append(df_sum)
# pd.concat(frames, sort=False).to_csv('/Users/oskar.hickl/binny_bench/img_hq_nc_bin_data_set_metrics_summary_v001.tsv', sep='\t')
#
# raise Exception
########################################################################################################################

for binner in set(df2.Binner.to_list()) | set('A'):
    df2 = df2.append(pd.Series([binner + '_placeholder' if col == 'Binner' else np.NAN for col in sorted(plot_data_dict.keys())],
                         index=sorted(plot_data_dict.keys()), name='Placeholder'))

binner_order = sorted(list(set(df2['Binner'].to_list())))

PROPS = {'boxprops': {'edgecolor': 'black'},
         'medianprops': {'color': 'black'},
         'whiskerprops': {'color': 'black'},
         'capprops': {'color': 'black'}}

qual_colors = ['#FFFFFF' if 'placeholder' in binner else 'whitesmoke' if 'HQ' in binner else 'silver' for
                             binner in binner_order]

ax = sns.violinplot(x="Recall", y="Binner", data=df2, width=1.9, fliersize=0,
                 linewidth=0.75, order=binner_order, palette=qual_colors,
                    scale="count", inner="quartile", **PROPS)

ax.tick_params(axis='y', pad=30)

plt.rcParams['legend.title_fontsize'] = 5

legend_palette = ['whitesmoke', 'silver']
fake_list = []
for e in legend_palette:
    fake_dot = ax.scatter(-10, 0, s=50, color=e, marker='o', linewidths=0.3, edgecolors='black')
    fake_list.append(fake_dot)
leg1 = ax.legend(fake_list, ['HQ', 'NC'], title='Bin Quality', facecolor='white', loc='upper right',
                 edgecolor='black', fontsize=4.5, markerscale=0.50)

lines = ax.get_lines()
categories = ax.get_yticks()[1::2]

median_x_pos = -0.1
if data_set_used == img_checkm_out_tsv:
    bin_count_x_pos = -0.05
else:
    bin_count_x_pos = -0.05

ax.text(bin_count_x_pos - 0, -0.6, 'Avg Recall [bp]', ha='center', va='center', fontweight='bold', size=5, color='black')

label_binner_order = [binner for binner in binner_order if 'placeholder' not in binner]

for cat, binner in zip(categories, label_binner_order):
    binner_avg_recall_sample_data = [str(i) for i in stat_dict_used[binner.split('_')[0]][binner.split('_')[1]]['Average Recall per Sample [bp]']]
    print(cat, binner, binner_avg_recall_sample_data)
    ax.text(bin_count_x_pos - 0, cat, binner_avg_recall_sample_data[0] + u"\u2009\u00B1\u2009" + binner_avg_recall_sample_data[1] + '%', ha='center',
            va='center', fontweight='semibold', size=5, color='black')

if data_set_used == img_checkm_out_tsv:
    data_set_palette = ['#45ac95', '#597bdb', '#6b7bbd', '#4cacd6', '#cf622f', '#8b5fd8', '#a15cba', '#b581c3', '#d352ae',
                        '#b9a333', '#bb5f92', '#9c8748', '#60b142', '#5b9f5d', '#d54156', '#c65b78', '#c4705a']
else:
    data_set_palette = ["#c9596f", "#709d60", "#7769d6", "#b68f35", "#b963b7", "#6895c8"]

data_set_hue_order = sorted([i for i in set(df2['Data Set'].to_list()) if type(i) is str])
if data_set_used == img_checkm_out_tsv:
    strip_jitter = 0.5
    x_limit_min = -1.25
else:
    strip_jitter = 0.5
    x_limit_min = -0.75

sns.stripplot(x="Recall", y="Binner", data=df2, size=3.2, color='white', linewidth=0.5, hue="Data Set",
              edgecolor='black', palette=data_set_palette, order=binner_order, dodge=True, jitter=strip_jitter, hue_order=data_set_hue_order)

ax.set_xlim(-0.01, 1)
plt.xticks(fontsize=5.5)

tick_spacing = 0.05
ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

label_list = []
for t in ax.get_legend_handles_labels():
    label_list.append(t)

new_list = []
for txt in label_list[1]:
    new_list.append(txt)
label_list[1] = new_list

if data_set_used == img_checkm_out_tsv:
    labels = ['Binny', 'CONCOCT', 'Maxbin2', 'MetaBAT2', 'VAMB']
    ax.set_yticks([pos + 1.0 for pos in range(len(binner_order) - 1)[1::4]])
    ax.set_yticklabels(labels, fontsize=7)
else:
    labels = ['Binny', 'CONCOCT', 'Maxbin2', 'MetaBAT2', 'VAMB']
    ax.set_yticks([pos + 1.0 for pos in range(len(binner_order) - 1)[1::4]])
    ax.set_yticklabels(labels, fontsize=7)

for vals in [pos + 0.0 for pos in range(len(binner_order) - 1)[4::4]]:
    y_ax_line = mlines.Line2D([-50, 250], [vals, vals], color='white', linewidth=1)
    ax.add_line(y_ax_line)
    y_ax_line2 = mlines.Line2D([-50, 250], [vals - 2, vals - 2], color='white', linewidth=0.25)
    ax.add_line(y_ax_line2)

if data_set_used == img_checkm_out_tsv:
    leg = plt.legend(handles=label_list[0],
               labels=label_list[1],
               fontsize=3.8, title='Environment', ncol=1, markerscale=0.50, handletextpad=0.5,
               labelspacing=1.15, columnspacing=1.5, edgecolor='black', facecolor='white')
else:
    leg = plt.legend(handles=label_list[0],
               labels=[' '.join([label_part.split('_')[0][:-1].upper(), label_part.split('_')[0][-1], label_part.split('_')[-1].title()]) for label_part in label_list[1]],
               fontsize="x-small", title='Data Set', edgecolor='black', facecolor='white', title_fontsize="small")
    print(label_list[1])

ax.add_artist(leg1)

plt.draw()

p = leg.get_window_extent()

if data_set_used == img_checkm_out_tsv:
    out_file_prefix = 'img'
    plt.savefig(
        '/Users/oskar.hickl/binny_bench/amber_data/benchmark_plots/{0}_genomes_recall_min{1}_{2}.pdf'.format(
            out_file_prefix, min_size, '_'.join([data_set_to_plot.split('\n')[0], data_set_to_plot.split('\n')[1]]).replace(':', '')),
        bbox_inches='tight')
else:
    out_file_prefix = 'cami'
    plt.savefig(
        '/Users/oskar.hickl/binny_bench/amber_data/benchmark_plots/{0}_genomes_recall_min{1}_max{2}_{3}.pdf'.format(
            out_file_prefix, min_size, max_size, data_set_to_plot),
        bbox_inches='tight')
plt.close()
