import re
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from functions import params
from functions import data_processing as dp
from functions import parsing

CER_GENOME = params.CER_GENOME
LAB_WT_NORM = params.LAB_WT_NORM
GP = params.GP

def find_motif_in_seq(seq, motif_regexp):
    rc_seq = str(Seq(seq).reverse_complement())
    mot_loc = []
    if re.search(motif_regexp, seq) != None:
        [mot_loc.append(i.start()+round(len(i.group())/2)) for i in re.finditer(motif_regexp, seq)]
    if re.search(motif_regexp, rc_seq) != None:
        [mot_loc.append(len(seq)-i.start()-1-round(len(i.group())/2)) for i in re.finditer(motif_regexp, rc_seq)]
    return mot_loc


def get_prom_seq(chromosome_str, start_pos, end_pos, prom_length):
    if start_pos < end_pos:
        seq = CER_GENOME[chromosome_str][start_pos-prom_length:start_pos]
    else:
        seq = str(Seq(CER_GENOME[chromosome_str][start_pos:start_pos+prom_length]).reverse_complement())
    return seq


def get_prom_signal(tf, chromosome, start_pos, end_pos, prom_length, motif_wind):
    if start_pos < end_pos:
        signal = LAB_WT_NORM[tf][chromosome-1][start_pos-prom_length-motif_wind:start_pos+motif_wind]
    else:
        signal = LAB_WT_NORM[tf][chromosome-1][start_pos+prom_length+motif_wind:start_pos-motif_wind:-1]
    return signal


def tf_enrichment_rank(pos, tested_tf, chromosome, start_pos,end_pos, prom_length, motif_wind, tf_list):
    tfs_wind_prc_signal = []
    for curr_tf in tf_list:
        tf_signal = np.array(get_prom_signal(curr_tf, chromosome, start_pos, end_pos, prom_length, motif_wind))
        tfs_wind_prc_signal.append(np.sum(tf_signal[pos-motif_wind:pos+motif_wind]) / np.sum(tf_signal))
    ordered_tfs = tf_list[np.argsort(tfs_wind_prc_signal)[::-1]]
    tested_tf_i = np.where(ordered_tfs==tested_tf)[0]
    return tested_tf_i


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def comb_proms(tf1, tf2, sum_prom, threshold):
    '''This function returns the union of bound promoters of a pair of tfs'''
    tf1_ind = sum_prom.loc[sum_prom[tf1]>threshold,tf1].sort_values(ascending=False).index
    tf2_ind = sum_prom.loc[sum_prom[tf2]>threshold,tf2].sort_values(ascending=False).index
    ids = pd.Index(np.concatenate((tf1_ind,tf2_ind))).drop_duplicates()
    return ids


def get_signal(curr_prom, tf, pad_signal, prom_length):
    chromosome = int(GP.iloc[curr_prom, :]['Chromosome'])
    start_pos = int(GP.iloc[curr_prom, :]['TSS_stein_Start'])
    end_pos = int(GP.iloc[curr_prom, :]['TSS_stein_End'])
    if start_pos < end_pos:
        signal = LAB_WT_NORM[tf][chromosome-1][start_pos-prom_length-pad_signal:start_pos+pad_signal]
    else:
        signal = LAB_WT_NORM[tf][chromosome-1][start_pos+prom_length+pad_signal:start_pos-pad_signal:-1]
    return signal


def find_motif_loc(seq, motif_regexp):
    rc_seq = str(Seq(seq).reverse_complement())
    mot_loc_f = []
    mot_loc_r = []
    if re.search(motif_regexp, seq) != None:
        [mot_loc_f.append(i.start()+round(len(i.group())/2)) for i in re.finditer(motif_regexp, seq)]
    if re.search(motif_regexp, rc_seq) != None:
        [mot_loc_r.append(len(seq)-i.start()-1-round(len(i.group())/2)) for i in re.finditer(motif_regexp, rc_seq)]
    l = mot_loc_f+mot_loc_r
    if sum(isinstance(i, list) for i in l) > 0:
        print(l)
        locs = [item for sublist in l for item in sublist]
    else:
        locs = l
    return locs


def get_bound_mot_locs(mot_locs,tf_signal,pad,window,quantile, prom_length):
    mot_per_sig = []
    for loc in mot_locs:
        mot_per_sig.append(np.sum(tf_signal[pad+loc-window:loc+pad+window]))
    quan_val = get_quantile_val(tf_signal, quantile, window, prom_length)
    bound_locs = list(np.array(mot_locs)[np.array(mot_per_sig)>quan_val])
    return bound_locs, mot_per_sig


def get_quantile_val(tf_signal, quantile, window, prom_length):
    wind_sum = []
    jumps = window*2
    for i in range(0, prom_length, jumps):
        wind_sum.append(np.sum(tf_signal[i:i+(window*2)-1]))
    quan_val = np.quantile(wind_sum, quantile)
    return quan_val


# def get_seq(curr_prom, prom_length):
#         chromosome = int(GP.iloc[curr_prom, :]['Chromosome'])
#         start_pos = int(GP.iloc[curr_prom, :]['TSS_stein_Start'])
#         end_pos = int(GP.iloc[curr_prom, :]['TSS_stein_End'])
#         chromosome_str = 'chr'+str(chromosome)
#         if start_pos < end_pos:
#             seq = CER_GENOME[chromosome_str][start_pos-prom_length:start_pos]
#         else:
#             seq = str(Seq(CER_GENOME[chromosome_str][start_pos:start_pos+prom_length]).reverse_complement())
#         return seq

def get_seq(curr_prom,prom_length):
    chromosome = int(GP.iloc[curr_prom,:]['Chromosome'])
    chromosome_str = 'chr'+str(chromosome)
    if np.isnan(GP.iloc[curr_prom,:]['TSS_stein_Start']):
        start_pos = int(GP.iloc[curr_prom,:]['ORF_Start'])
        end_pos = int(GP.iloc[curr_prom,:]['ORF_End'])
    else:
        start_pos = int(GP.iloc[curr_prom,:]['TSS_stein_Start'])
        end_pos = int(GP.iloc[curr_prom,:]['TSS_stein_End'])
    if start_pos < end_pos:
        seq = CER_GENOME[chromosome_str][start_pos-prom_length:start_pos]
    else:
        seq = str(Seq(CER_GENOME[chromosome_str][start_pos:start_pos+prom_length]).reverse_complement())
    return seq