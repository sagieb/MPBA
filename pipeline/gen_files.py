import pandas as pd
import sys
import os
import re

def get_sample_info(file_path,samp_info,plasmid_info, output):
    samp_info = pd.read_csv(samp_info)
    plasmid_info = pd.read_csv(plasmid_info)   
    well_barcode = int(file_path.split('.')[1])
    pool_name = file_path.split('.')[0].split('/')[1].split('_')
    pool_name = '_'.join(pool_name[:3])
    tf = samp_info.loc[(samp_info['Pool name']==pool_name) & (samp_info['Well barcode number']==well_barcode)]['Strain'].values[0]
    pool_number = samp_info.loc[(samp_info['Pool name']==pool_name) & (samp_info['Well barcode number']==well_barcode)]['Plasmid pool number'].values[0]
    prom_nums = list(plasmid_info.loc[(plasmid_info['Factor']==tf) & (plasmid_info['Pool']==pool_number)]['Plasmid'])
    prom_names = ['prom_' + str(name) for name in prom_nums]
#     prom_names.append('orf_abs')
#     print(prom_names)
    with open(output, 'a') as f:
        for prom in prom_names:
            f.write(file_path+','+prom+'\n')
    return prom_names


if __name__ == "__main__":
    prom_names = get_sample_info(sys.argv[1],sys.argv[2],sys.argv[3], sys.argv[4])
#     for prom in prom_names:
#         mutation_count(sys.argv[1], prom)
#     with open(sys.argv[4], 'w') as f:
#         f.write(' '.join(prom_names))