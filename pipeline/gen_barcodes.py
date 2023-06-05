import pandas as pd
import yaml
import io
import os
import glob
import re
from pathlib import Path

if __name__ == "__main__":
    
    g1 = []
    g2 = []
    g3 = []
    full = []
    for name in os.listdir('raw_data/'):
        r = re.search("(.*)?_S(.*?)_R(.*?)\.fastq\.gz", name)
        g1.append(r.group(1))
        g2.append(r.group(1)+'_S'+r.group(2))
        # g3.append(r.group(3))
        # full.append(name)

    connected = dict(zip(g1,g2))

    tab = pd.read_csv('libs_info/samp_info.csv', index_col=0)
    barcodes = dict()
    for exp in connected.keys():
        print(exp)
        barcodes[connected[exp]] = list(map(str, list(tab.loc[tab.loc[:, 'Pool name'] == exp].loc[:, 'Well barcode number'])))

    yaml_config = {
    'barcodes': 'barcodes/barcodes.txt',
    'prom_adapt':'libs_info/prom_adapt.fasta',
    'samp_info': 'libs_info/samp_info.csv',
    'plasmid_info': 'libs_info/plasmids_by_pools.csv',
    'adaptor': 'libs_info/bcs_cut.fasta',
    'sample':'',
    'personal_barcodes':'',
    'end_pos_check':'131',
    'length_to_check':'6',
    'promoter_sequences':'libs_info/promoter_seqs.csv'
    }

    yaml_config['personal_barcodes'] =barcodes 

    for name in connected.values():
        per_conf = yaml_config
        per_conf['sample'] = name
        path = Path(name+'.yaml')
        if path.is_file():
            print('exists')
        else:
            with io.open('{}.yaml'.format(name), 'w') as outfile:
                yaml.dump(per_conf, outfile, default_flow_style=False, allow_unicode=True)