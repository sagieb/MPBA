import numpy as np
import pandas as pd
import sys
import os
from functions import parsing
from functions import data_processing as dp
import scipy.io as sio

# used paths
curr_dir = os.getcwd()
data_path = os.path.join(curr_dir, 'data')  # directory containing all data used
genome_annot_path = os.path.join(data_path, 'genome_annot')  # a directory containing all genomic information
libs_info_path = os.path.join(data_path, 'libs_info')  # directory containing libraries information
lab_binding_path = os.path.join(data_path, 'lab_binding_data')
lib_data_path = os.path.join(data_path, 'libs_data')

# Library data and results
ALL_TF_COOP = pd.read_csv(os.path.join(lib_data_path, 'all_tf_coop.csv'))  # all pairwise cooperativity
# scores (including TF motif)
RES_PATH = os.path.join(lib_data_path, 'res_files')  # directory containing all filtered resuslt tables in csv files

# Library information files and annotation folder
LIBS_ANNOT_PATH = os.path.join(libs_info_path, 'libs_annotations')  # containing all libraries annotations files
LIBS_RANK = pd.read_csv(os.path.join(libs_info_path, 'rank_libs.csv'))
LIBS_RATIO = pd.read_csv(os.path.join(libs_info_path, 'ratiomax_libs.csv'))
LIBS_INFO = pd.read_csv(os.path.join(libs_info_path, 'all_libs_TS.csv'))
LIBS_RATIO.dropna(how='all', inplace=True)

# genomic information files
GP = pd.read_csv(os.path.join(genome_annot_path, 'gene_info.csv'), index_col=0)
CER_GENOME = parsing.decompress_pickle(os.path.join(genome_annot_path, 'cerGenome.pbz2'))
PROM_POS = np.load(os.path.join(genome_annot_path, 'promoter_pos_1000_150_into.npy'))
gene_name = sio.loadmat(os.path.join(genome_annot_path, 'gene_names.mat'), simplify_cells=True)['gene_names']
GENE_NAMES = gene_name.tolist()
STEIN = np.load(os.path.join(genome_annot_path, 'stein.npy'))
PROMOTER_LENGTH_GP = np.load(os.path.join(genome_annot_path, 'prompter_length_GP.npy'))

# genomic binding signal files
LAB_WT_NORM = parsing.load_dict(os.path.join(lab_binding_path, 'norm_chec_by_factor'))
LAB_WT_ZSCORE = parsing.load_dict(os.path.join(lab_binding_path, 'zcore_chec_by_factor'))
TF_NAMES = list(LAB_WT_NORM.keys())
TF_MOTIFS = pd.read_csv(os.path.join(lab_binding_path, 'tf_motifs.csv'), index_col=0)
TF_SUMPROM = pd.read_csv(os.path.join(lab_binding_path, 'tfs_sumprom.csv'))
MOT_EUC_DIST = pd.read_csv(os.path.join(lab_binding_path, 'euclidean_dist.csv'), index_col=0)
MNASE = sio.loadmat(os.path.join(lab_binding_path, 'MNase_OD4.mat'))['MNase_OD4'][0]

# Fixed parameters
PROMOTER_LENGTH = 1000
TOP_PROM_NUM = 100
LIB_SIZE = 164
INTO_GENE = 150
PRIMER_LENGTH = 18


if __name__ == '__main__':
    x = 1
