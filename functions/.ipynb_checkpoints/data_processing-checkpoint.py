import os
import numpy as np
import pandas as pd
from glob import glob
import re
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
import scipy.io
import warnings
from matplotlib import gridspec
from matplotlib.cm import ScalarMappable
from pptx import Presentation
from pptx.util import Inches
import io
from mpl_toolkits.mplot3d import Axes3D
sys.path.insert(0,'/home/labs/barkailab/tamarj/TF_combinatorics/')
from functions import parsing
from functions import promoter_analysis

data_path = '/home/labs/barkailab/tamarj/TF_combinatorics/data' # directory containing all project data and info
RES_PATH = os.path.join(data_path, 'method_exps_res') # directory containing all result tables in csv files
libs_info_path = os.path.join(data_path, 'libs_info') # directory containing libraries informaiton

LIBS_RANK = pd.read_csv(os.path.join(libs_info_path, 'rank_libs.csv'))
LIBS_RATIO = pd.read_csv(os.path.join(libs_info_path, 'ratiomax_libs.csv'))
LIBS_INFO = pd.read_csv(os.path.join(libs_info_path, 'all_libs_TS.csv'))
LIBS_RATIO.dropna(how='all',inplace=True)

CER_GENOME = parsing.decompress_pickle('/home/labs/barkailab/tamarj/pythonProject/LabUtils/lab_utils/data/cerGenome/cerGenome.pbz2')
# LAB_WT_NORM = parsing.mat2py('/home/labs/barkailab/tamarj/CB/data/lab_wt_norm.mat')
PROM_POS = np.load('/home/labs/barkailab/tamarj/pythonProject/LabUtils/lab_utils/data/promoter_pos_1000_150_into.npy')
GENE_NAMES = list(promoter_analysis.GENE_NAMES)

def get_lib_info(lib_num):
    ## this function gives useful information about the library. output is a dictionary.
    mut_dict = {'R': '[GA]', 'M': '[AC]', 'K': '[GT]'}  # defining a dict for mapping our variable position
    info = {} # building and info dict containing useful information about the current library
    curr_lib_info = lib_info.loc[lib_info["Number"] == int(lib_num)]  #library numnber
    info['gene'] = curr_lib_info["Gene"].values[0]  # the library is sequence was taken from this gene's promoter
    info['gene_id'] = gene_names.index(gene.upper())  # getting the gene ID (based on GP)
    info['lib_seq'] = curr_lib_info['Library_sequence'].values[0]  #library sequnce (with variable position nuclotides)
    info['wt_seq'] = curr_lib_info['Sequence'].values[0]  #Genomic sequence the library was based on
    info['mut_by_tf'] = curr_lib_info['TFs_by_motif'].values[0].split(',')  # List of tfs for which we mutated the motifs
    info['start_position'] = int(curr_lib_info["Promoter_start_position"].values[0])  #position in the promoter where the library is located
    info['mut_loc'] = [int(m.start()) for m in re.finditer('[RYMKSWHBVDN]', info['lib_seq'])]  # position within the library sequence where mutations are located
    info['wt_at_var_loc'] = [info['wt_seq'][loc] for loc in info['mut_loc']]  # wild type/genomic nucleotides found in variable positions
    info['var_pos_pattern'] = ''.join([mut_dict[info['lib_seq'][i]] for i in info['mut_loc']])  #string containing both options wt/mut in each var position
    return info

    def filter_results(tf, lib_num, exp):
        # this function generates a table with the number of reads recieved for true sequences
        # (containing a combination of true ordered nucleotides in variable positions)
        curr_lib_info = get_lib_info(lib_num)  # using get_lib_info to get needed information
        file_name = tf + '_' + lib_num + '_' + exp + '.csv'  # looking for the relevant results file based on our system ('tf_libnum_exp.csv')
        res_table = (pd.read_csv(os.path.join(RES_PATH, file_name), index_col=0))  # reading current result table into pd dataframe
        res_table.index.name = 'var_sequence'  # naming the index column (which is currently empty)
        res_table_true_seqs = res_table.iloc[np.where(res_table.index.str.findall(curr_lib_info['var_pos_pattern']))[0], :]
        # subseting the res table taking only location with true sequences
        return res_table_true_seqs