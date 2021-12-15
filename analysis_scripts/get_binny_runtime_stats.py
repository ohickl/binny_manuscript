import glob
import pandas as pd

sample_path_list = glob.glob('binny_out/c?_t_*/*')

runtime_stats_list = []
step_dict = {'annotate': 'Prokka', 'checkm': 'Mantis', 'binny': 'binny'}
df_header = ['data_set', 'sample', 'step', 'h:m:s', 'm']

for sample_path in sample_path_list:
    data_set = sample_path.split('/')[1]
    sample = sample_path.split('/')[2]
    sample_benchmark_files = glob.glob(sample_path + '/logs_original/*_benchmark.txt')
    if sample_benchmark_files:
        total_time = 0
        for bench_file in sample_benchmark_files:
            step = bench_file.split('/')[-1].split('_')[1]
            step = step_dict[step]
            with open(bench_file, 'r') as f:
                next(f)
                for line in f:
                    time_taken = line.split('\t')[1]
                    time_taken_min = round(float(line.split('\t')[0]) / 60, 1)
                    total_time += time_taken_min
                    # print('\t'.join([data_set, sample, step, time_taken, str(time_taken_min)]))
            runtime_stats_list.append([data_set, sample, step, time_taken, time_taken_min])
        runtime_stats_list.append([data_set, sample, 'total', 'NaN', total_time])

binny_runtime_stats = pd.DataFrame(runtime_stats_list, columns=df_header)
binny_runtime_stats = binny_runtime_stats.sort_values(by=['data_set', 'sample', 'step'])

binny_runtime_stats.to_csv('binny_out/binny_runtime_stats.tsv', sep='\t')

for i in ['total', 'Prokka', 'Mantis', 'binny']:
    step_binny_runtime_stats = binny_runtime_stats[binny_runtime_stats['step'] == i].describe(include='all')
    print(step_binny_runtime_stats)
    step_binny_runtime_stats.to_csv('binny_out/binny_runtime_stats_summary_{0}.tsv'.format(i), sep='\t')

