import os
import pandas as pd
import numpy as np
import re
import itertools

'''This script generates the final csv files containing the data for all measured time points and repeats for a given  
TF in a given experiment.'''


def get_all_combinations(pattern_str):
    char_combinations = re.findall(r'\[([^]]+)\]', pattern_str)  # Find all character combinations within square brackets
    combinations = list(itertools.product(*char_combinations))  # Generate all combinations
    combinations = [''.join(combo) for combo in combinations]  # Convert combinations to strings
    return combinations

# generating relevant paths and loading data


run_path = './'
MBPA_path = '../'
res_folder = 'final_files'
general_data_folder = 'data'
lib_info_folder = 'libs_info'
sample_info_file = 'samp_info.csv'
libs_info_file = 'all_libs_TS.csv'
libs_data_dir = 'libs_data'
csv_dir = 'res_files'

res_data_path = os.path.join(run_path, res_folder)
general_data_path = os.path.join(MBPA_path, general_data_folder)
libs_info_path = os.path.join(general_data_path, lib_info_folder)
libs_data_path = os.path.join(general_data_path, libs_data_dir)

sample_info = pd.read_csv(os.path.join(libs_info_path, sample_info_file))
libs_info = pd.read_csv(os.path.join(libs_info_path, libs_info_file), index_col=[0])
data_files = os.listdir(res_data_path)

mut_dict = {'R': '[GA]', 'M': '[CA]', 'K': '[GT]', 'Y': '[ST]'}
sample_info.columns = [col.replace(' ', '_') for col in sample_info.columns]
relevant_files = [file for file in data_files if '.txt' in file]

# generate the directory in which the final csv files will be saved.
out_path = os.path.join(libs_data_path, csv_dir)
os.mkdir(out_path)

# Gather the relevant results for each TF in all time points and repeats for a given experiment:
for file in relevant_files:
    if file not in relevant_files:
        continue
    curr_prom_number = file.split('.')[-2].split('_')[-1]
    sample_name = file.split('.')[0]
    pool_name = sample_name.split('_S')[0]
    well_barcode = int(file.split('.')[1])
    curr_samp_info = sample_info.query("Pool_name == @pool_name & Well_barcode_number == @well_barcode")
    curr_exp = curr_samp_info['Exp'].values[0]
    curr_tf = curr_samp_info['Strain'].values[0]
    curr_plasmid_pool = curr_samp_info['Plasmid_pool_number'].values[0]
    res_subset = sample_info.query("Strain == @curr_tf & Plasmid_pool_number == @curr_plasmid_pool & Exp == @curr_exp")
    res_files = []
    for row_i in res_subset.index:
        curr_row = res_subset.loc[row_i]
        composed_name = re.compile(curr_row['Pool_name'] + '_S[0-9]{1,2}' + '.' + str(curr_row['Well_barcode_number']) +
                                   '.prom_' + curr_prom_number + '.txt')
        res_files.append(list(filter(composed_name.match, relevant_files))[0])
    # generate all possible library sequences:
    curr_lib_info = libs_info.loc[int(curr_prom_number)]
    curr_lib_seq = curr_lib_info['Library_sequence']
    pattern = []
    for char in curr_lib_seq:
        if char in mut_dict.keys():
            pattern.append(mut_dict[char])
    str_pattern = ''.join(pattern)
    curr_lib_seqs = get_all_combinations(str_pattern)
    curr_res_table = pd.DataFrame(index=curr_lib_seqs)
    for file_i, curr_file in enumerate(res_files):
        curr_info = res_subset.iloc[file_i]
        curr_col_name = '_'. join([str(curr_info['Trans_repeat']), str(curr_info['Time_point']),
                                   str(curr_info['Time_point_repeat'])])
        curr_file_data = pd.read_csv(os.path.join(res_data_path, curr_file), sep=r'\s+', header=None,
                                     names=['read_num', 'seq'])
        read_counts = []
        for seq in curr_lib_seqs:
            if seq in curr_file_data['seq'].values:
                read_counts.append(curr_file_data.loc[curr_file_data['seq'] == seq]['read_num'].values[0])
            else:
                read_counts.append(0)
        curr_res_table[curr_col_name] = read_counts
    curr_res_table_ordered = curr_res_table.sort_index(axis=1)
    csv_name = '_'.join([curr_tf, curr_prom_number, str(curr_exp)])
    curr_res_table_ordered.to_csv(os.path.join(out_path, csv_name+'.csv'))

    relevant_files = np.setdiff1d(relevant_files, res_files)
