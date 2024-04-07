import os
import numpy as np
import pandas as pd
import re
import scipy.io as sio
import sys
from functions import parsing
from functions import prom_analysis
from functions import params
from glob import glob
import pickle
from Bio import SeqIO


def get_lib_info(lib_num):
    ''' This function gives useful information about the library. output is a dictionary.'''
    mut_dict = {'R': '[GA]', 'M': '[AC]', 'K': '[GT]'}  # defining a dict for mapping our variable position
    info = {} # building and info dict containing useful information about the current library
    curr_lib_info = params.LIBS_INFO.loc[params.LIBS_INFO["Number"] == int(lib_num)]  #library numnber
    info['gene'] = curr_lib_info["Gene"].values[0]  # the library is sequence was taken from this gene's promoter
    info['gene_id'] = params.GENE_NAMES.index(info['gene'].upper())  # getting the gene ID (based on GP)
    info['lib_seq'] = curr_lib_info['Library_sequence'].values[0]  #library sequnce (with variable position nuclotides)
    info['wt_seq'] = curr_lib_info['Sequence'].values[0]  #Genomic sequence the library was based on
    info['mut_by_tf'] = curr_lib_info['TFs_by_motif'].values[0].split(',')  # List of tfs for which we mutated the motifs
    info['start_position'] = int(curr_lib_info["Promoter_start_position"].values[0])  #position in the promoter where the library is located
    info['mut_loc'] = [int(m.start()) for m in re.finditer('[RYMKSWHBVDN]', info['lib_seq'])]  # position within the library sequence where mutations are located
    info['wt_at_var_loc'] = [info['wt_seq'][loc] for loc in info['mut_loc']]  # wild type/genomic nucleotides found in variable positions
    info['var_pos_pattern'] = ''.join([mut_dict[info['lib_seq'][i]] for i in info['mut_loc']])  #string containing both options wt/mut in each var position
    return info


def sort_table_by_tp(res_table):
    '''This function orders a given resulte table (seqs x samples) by tps, according to columns names.'''
    ord_ids = np.argsort([int(x.split('_')[1]) for x in list(res_table.columns)])
    ord_table = res_table.iloc[:, ord_ids]
    return ord_table


def filter_results(tf, lib_num, exp):
    ''' This function generates a table with the number of reads recieved for true sequences
    (containing a combination of true ordered nucleotides in variable positions).'''
    curr_lib_info = get_lib_info(lib_num)  # using get_lib_info to get needed information
    file_name = tf + '_' + str(lib_num) + '_' + str(exp) + '.csv'  # looking for the relevant results file based on our system ('tf_libnum_exp.csv')
    res_table = (pd.read_csv(os.path.join(params.RES_PATH, file_name), index_col=0))  # reading current result table into pd dataframe
    res_table.index.name = 'var_sequence'  # naming the index column (which is currently empty)
    res_table_true_seqs = res_table.iloc[np.where(res_table.index.str.findall(curr_lib_info['var_pos_pattern']))[0], :]  # subseting the res table taking only location with true sequences
    reordered_filterd_table = sort_table_by_tp(res_table_true_seqs)
    return reordered_filterd_table

def norm_reads(filtered_results):
    '''This function gets a table with the number of reads recieved for each *true* sequence and normalizes to
     "norm_to". Then, bad samples are filtered out using the rm_samples function.'''
    norm_to = 10**6
    normalized_results = filtered_results/filtered_results.sum(axis=0)*norm_to
    return normalized_results


def res_log2(normalized_results):
    '''this function converts your normalized results to log2.'''
    log2_results = (normalized_results+1).apply(np.log2)
    return log2_results


def get_samp_info(res_table):
    ''' This function takes the column names of a results table and generates a matrix where
    each row is a sample, col1 = biological rep, col2 = time-point, col3 = technical rep.'''
    samp_inf = [col_name.split('_') for col_name in list(res_table)]
    sample_info_mat = np.asarray(samp_inf,np.ndarray).astype(int)
    return sample_info_mat


def mean_bio_rep_tps(norm_res, sample_info):
    '''This function means over tps within each biological repeat.'''
    bio_repeat_id = 0
    tp_id = 1  # indices by sample info: col1 = biological rep, col2 = time-point, col3 = technical rep.
    bio_repeats = list(np.unique(sample_info[:, bio_repeat_id]))  # get biological repeats and biological repeats number
    num_bio_repeats = len(bio_repeats)
    tps = list(np.unique(sample_info[:, tp_id]))  #get time points and time points number
    num_tps = len(tps)
    samples_tps = np.sort(tps*num_bio_repeats)  # list of time points for each biological repeat
    samples_bio_repeats = bio_repeats*num_tps  # list of biological repeats for each time point
    mean_col_names = ['_'.join([str(curr_bio), str(curr_tp)])+'_' for curr_bio, curr_tp in zip(samples_bio_repeats, samples_tps)]  # generate columns names for new mean datafram consisting of the relevant biological repeat and tp
    mean_repeats = pd.DataFrame()  # empty dataframe for mean biological repeats tps
    for i, col in enumerate(mean_col_names):  # iterate mean df columns names that
        curr_mean = norm_res.iloc[:, norm_res.columns.str.contains(col)].mean(axis=1)  # from norm_res mean over columns that therir column name contains the relevant biological repeat and tp
        mean_repeats = pd.concat([mean_repeats, curr_mean.rename(mean_col_names[i]).to_frame()],axis=1)  # add the mean data to the new mean dataframe
    mean_repeats.columns = [col[:-1] for col in mean_repeats.columns]
    return mean_repeats


def sem_bio_rep_tps(norm_res, sample_info):
    '''This function means over tps within each biological repeat.'''
    bio_repeat_id = 0
    tp_id = 1  # indices by sample info: col1 = biological rep, col2 = time-point, col3 = technical rep.
    bio_repeats = list(np.unique(sample_info[:, bio_repeat_id]))  # get biological repeats and biological repeats number
    num_bio_repeats = len(bio_repeats)
    tps = list(np.unique(sample_info[:, tp_id]))  #get time points and time points number
    num_tps = len(tps)
    samples_tps = np.sort(tps*num_bio_repeats)  # list of time points for each biological repeat
    samples_bio_repeats = bio_repeats*num_tps  # list of biological repeats for each time point
    sem_col_names = ['_'.join([str(curr_bio), str(curr_tp)])+'_' for curr_bio, curr_tp in zip(samples_bio_repeats, samples_tps)]  # generate columns names for new sem datafram consisting of the relevant biological repeat and tp
    sem_repeats = pd.DataFrame()  # empty dataframe for sem biological repeats tps
    for i, col in enumerate(sem_col_names):  # iterate sem df columns names that
        curr_sem = norm_res.iloc[:, norm_res.columns.str.contains(col)].sem(axis=1)  # from norm_res sem over columns that therir column name contains the relevant biological repeat and tp
        sem_repeats = pd.concat([sem_repeats, curr_sem.rename(sem_col_names[i]).to_frame()],axis=1)  # add the sem data to the new sem dataframe
    sem_repeats.columns = [col[:-1] for col in sem_repeats.columns]
    return sem_repeats


def norm_to_tp_0(mean_repeats_df, sample_info, space=None):
    ''' This function normalizes your activated samples to time zero
    if your data frame is on linear scale you should indicate 'linear' under the "space" input.'''
    mean_activate = mean_repeats_df.iloc[:, ~mean_repeats_df.columns.str.contains('_0')] # get sub table of only activated samples
    tps = list(np.unique(sample_info[:,1]))  #get time points and time points number
    num_tps = len(tps)
    mean_zero = mean_repeats_df.iloc[:,mean_repeats_df.columns.str.contains('_0')] # get sub table of only 0 tps
    replicate_mean_zero = pd.concat([mean_zero] * (num_tps-1), axis=1, ignore_index=True) # replicate time zero dataframe according to number of activated time points
    replicate_mean_zero.columns = mean_activate.columns # keep relevant column names
    if space == 'linear': #if your table is on linear space you divide in order to normalize to time 0. if it's on log2 scale you need to substract.
        tp_0_norm = mean_activate.div(replicate_mean_zero, axis=1)
    else:
        tp_0_norm = mean_activate.sub(replicate_mean_zero, axis=1)
    return tp_0_norm


def find_lib_files(lib):
    '''This function takes a library number and returns the names of all results files of this library.'''
    res_files = os.listdir(params.RES_PATH) # get a list of all result files
    lib_files = [file for file in res_files if '_'+str(lib)+'_' in file] # getting relevant library files out of all res files
    return lib_files


def find_tf_files(tf):
    '''This function takes a tf and returns the names of all results files of this tf.'''
    res_files = os.listdir(params.RES_PATH) # get a list of all result files
    lib_files = [file for file in res_files if str(tf)+'_' in file] # getting relevant tf files out of all res files
    return lib_files


def get_lib_tf_pool(tf, lib_num):
    '''This function gets TF name and library number and returns the library pool.'''
    pool_libs_info = pd.read_csv(os.path.join(params.libs_info_path, 'libs_by_pool.csv')) # load library by pools csv
    cond = (pool_libs_info["Factor"] == tf) & (pool_libs_info["Lib"] == int(lib_num)) # find the line matches the TF and library number
    rel_rows = pool_libs_info.loc[cond]
    if len(rel_rows) > 1: # in case we had only 1 pool of the top TF libraries the pool name is 'top', otherwise - returns the pool number
        return 'top'
    else:
        return rel_rows['Pool'].values[0]


def get_exp_absolute(tf, lib, exp):
    '''This function takes in a tf, library and experiment and returns a one line organized table with the number of absolute reads in each sample'''
    abs_table = params.ABS
    cond = ((abs_table["tf"] == tf) & (abs_table["promoter"] == int(lib)) & (abs_table["exp"] == int(exp))) # conditioning over the large
    # absolute table to get only lines that belong to the current exp,lib and tf
    sub_df = abs_table.loc[cond] #subset absolute table base on this condition
    col_names = ['%d_%d_%d' %(t.bio_rep, t.time_point, t.tech_rep) for _, t in sub_df.iterrows()] # generating column names to fit our general use
    abs_table_subset = pd.pivot_table(sub_df, values="abs_count", columns=col_names)  # pivot the table to get a line with all abs measurements
    abs_table_subset = sort_table_by_tp(abs_table_subset)  # sort the new table based on column names
    return abs_table_subset


def rm_samples(tf, lib_num, exp_num, res_table):
    '''This function removes from the result table the samples classified as bad. It gets the result table
    (before or after normalization) and the relevant tf, library number and experiment number.'''
    exp_num = int(exp_num)
    samples_to_remove = pd.read_csv(os.path.join(params.libs_info_path, 'samples_to_remove.csv'))
    # read csv containing info about samples to remove
    curr_samples_to_rm = samples_to_remove.loc[(samples_to_remove['Lib_number'] == int(lib_num)) &
                                               (samples_to_remove['Exp_number'] == exp_num) &
                                               (samples_to_remove['TF'] == tf)]
    # get the names of samples to remove relevant to the input
    max_tech_rep_num = max(
        [int(name.split('_')[-1]) for name in res_table.columns])  # number of time point repeats in the experiment
    var_num = len(res_table)  # number of variants
    samples_to_rm_names = curr_samples_to_rm['Sample_name']  # the sample names that will be removed

    sample_bio_rep = np.array(
        [int(name.split('_')[0]) for name in samples_to_rm_names])  # the bio repeats of the bad samples
    sample_tps = np.array(
        [int(name.split('_')[1]) for name in samples_to_rm_names])  # the time point of the bad samples

    bad_bio_rep = ''
    for i, samp in enumerate(samples_to_rm_names):  # iterate over the names of samples to remove
        curr_tp = sample_tps[i]  # current time point
        curr_bio_rep = sample_bio_rep[i]  # current bio repeat
        np.unique(res_table.columns.str.split('_')[0])
        curr_biorep_tp = str(curr_bio_rep) + '_' + str(curr_tp) + '_'
        same_tp_biorep_samples_rm = [s for i, s in enumerate(samples_to_rm_names.values) if
                                     curr_biorep_tp in samples_to_rm_names.values[i]]
        # get all sample names to remove with the current tp and bio repeat
        all_bio_rep_tp_samples = res_table.columns[res_table.columns.str.contains(curr_biorep_tp)]
        # get all sample names (good and bad samples) from same tp and bio repeat as current sample
        good_samples = [s for s in all_bio_rep_tp_samples if sum(samples_to_rm_names.str.contains(s)) == 0]
        # returns only samples from the same time point and bio repeat that are good

        if (len(same_tp_biorep_samples_rm) == 1) & (max_tech_rep_num == 2):
            # when removing only one tp repeat of the bio repeat and there are 2 time point repeats
            res_table[samp] = res_table[good_samples]  # change bad sample values to the other good time point repeat
        elif (len(same_tp_biorep_samples_rm) == 2) & (max_tech_rep_num == 2):
            # when removing 2 tp repeats of the bio repeat and there are 2 time point repeats
            res_table[samp] = np.nan * np.ones(shape=(1, var_num))[0]  # change bad sample values to nans
            bad_bio_rep = str(curr_bio_rep)  # in case one bio repeat has an entire tp to remove - save the bio repeat number
        elif (len(same_tp_biorep_samples_rm) == 1) & (max_tech_rep_num == 3):
            # when removing only one tp repeat of the bio repeat and there are 3 time point repeats
            res_table[samp] = res_table[res_table.columns.intersection(good_samples)].mean(axis=1)
            # change bad sample values to the mean of the other good time point repeats
        elif (len(same_tp_biorep_samples_rm) == 2) & (max_tech_rep_num == 3):
            # when removing 2 tp repeats of the bio repeat and there are 3 time point repeats
            res_table[samp] = res_table[good_samples]
            # change bad sample values to the one good time point repeat
        elif (len(same_tp_biorep_samples_rm) == 3) & (max_tech_rep_num == 3):
            # when removing 3 tp repeats of the bio repeat and there are 3 time point repeats
            res_table[samp] = np.nan * np.ones(shape=(1, var_num))[0]
            # change bad sample values to nans (no such scenario so far)
            bad_bio_rep = str(curr_bio_rep)
            # in case one bio repeat has an entire tp to remove - save the bio repeat number

    res_table_scatter = res_table.copy()  # copy the result table
    if bad_bio_rep:  # when one bio repeat has an entire tp to remove
        all_bio_rep = res_table_scatter.columns.str.split('_').str[0]  # get all the bio repeat sample names
        res_table_scatter[res_table_scatter.columns[all_bio_rep == bad_bio_rep]] = res_table_scatter[
            res_table_scatter.columns[all_bio_rep != bad_bio_rep]]
        # copy the "good" bio repeat instead of the "bad" ones
    return res_table, res_table_scatter


def mean_over_bio_reps(mean_tp_df):
    '''This function gets a data frame containing the mean over time point repeats and
    returns the mean over biological repeats'''
    col_names = mean_tp_df.columns # get column names
    cols_tps = col_names.str.split('_').str[1] # get time points
    tps = cols_tps.unique()
    mean_bio_rep = pd.DataFrame() # initiate a new data frame for storing the mean values
    for tp in tps: # iterate over time points
        mean_bio_rep[tp] = mean_tp_df[col_names[cols_tps == tp]].mean(axis=1)
    return mean_bio_rep


def sem_over_bio_reps(mean_tp_df):
    '''This function gets a data frame containing the mean over time point repeats and
    returns the mean over biological repeats'''
    col_names = mean_tp_df.columns # get column names
    cols_tps = col_names.str.split('_').str[1] # get time points
    tps = cols_tps.unique()
    sem_bio_rep = pd.DataFrame() # initiate a new data frame for storing the mean values
    for tp in tps: # iterate over time points
       sem_bio_rep[tp] = mean_tp_df[col_names[cols_tps == tp]].sem(axis=1)
    return sem_bio_rep


def norm_var_to_abs(tf,lib_num,exp):
    '''This function finds the absolute change of each sequence variant. technical repeats are averaged.'''
    file_name = tf + '_' + str(lib_num) + '_' + str(exp) + '.csv'  # find results file ('tf_libnum_exp.csv')
    filt_raw = (pd.read_csv(os.path.join(params.RES_PATH, file_name), index_col=0))
    filt_bad_samples, _ = rm_samples(tf,lib_num,exp, filt_raw)
    abs_raw, _ = rm_samples(tf, lib_num, exp, get_exp_absolute(tf,lib_num,exp)) # get absolute measurements of all samples
    sample_info = get_samp_info(filt_bad_samples) # get sample info to be used for averaging
    data_tech_mean = mean_bio_rep_tps(filt_bad_samples,sample_info) # average tech reps for each seq variant
    abs_tech_mean = mean_bio_rep_tps(abs_raw,sample_info) # average tech reps of absoute measurement
    abs_norm = data_tech_mean.div(abs_tech_mean.values).apply(np.log2) # normalize each variant to the absolute and transform to log2
    sample_info = get_samp_info(abs_norm) # get sample info to be used for tp0 normalization
    norm_abs_tp0 = norm_to_tp_0(abs_norm,sample_info) # normalize to timepoint 0
    return norm_abs_tp0


def get_norm_avg(tf,lib_num,exp):
    '''This function returns the averaged (tech & bio) data normalized to time0 after bad sample clean'''
    file_name = tf + '_' + str(lib_num) + '_' + str(exp) + '.csv'  # find results file ('tf_libnum_exp.csv')
    res_table = (pd.read_csv(os.path.join(params.RES_PATH, file_name), index_col=0))
    norm_data = dp.norm_reads(res_table)
    filt_norm, _ = rm_samples(tf, lib_num, exp, norm_data)
    sample_info = get_samp_info(filt_norm)
    norm_mean_tech = mean_bio_rep_tps(filt_norm, sample_info).apply(np.log2)
    sample_info = get_samp_info(norm_mean_tech)
    relative_tp0_norm = norm_to_tp_0(norm_mean_tech, sample_info)
    fin_data = mean_over_bio_reps(relative_tp0_norm)
    return fin_data


def get_motif_info_by_pos(lib):
    '''This function recieves a library number and generates a dataframe with TF
    motif annotation and the variable position within lib their motif contains'''
    if type(lib) != str: # if an integer is recieved convert to str
        lib = str(lib)
    mdir = params.LIBS_ANNOT_PATH # the directory containing library annotations from benchling
    lib_patt = '/Lib' + lib #pattern to use to locate a wanted annotation file (by lib)
    pf = 'CGATGCGCATGCGTACGC' # primter F to later remove length from start
    file = glob(mdir+lib_patt+'*') # find the wanted lib annotation
    with open(file[0]) as input_handle: # prase genebank annotation file
        record = list(SeqIO.to_dict(SeqIO.parse(input_handle, "genbank")).values())[0]
    tf_loc_cols = ['TF','Start','End'] # columns to be used in initial dataframe
    annot_df = pd.DataFrame() # initial df to contain all motif annotations
    for curr_feat in record.features: # run over features (motif annotations)
        if 'lift' not in curr_feat.qualifiers['label'][0]: # ignore primers called lift F/R
            labels = curr_feat.qualifiers['label'][0].split(' ')[0].split('/') # remove stuff after spaces and split found TFs on the same location seperated by /
            for i,label in enumerate(labels[1:]): # sometimes we wrote Fkh1/2. this is a fix. the list will be ['Fkh1','Fkh2']
                if len(label)==1:
                    labels[i+1] = labels[0][0:-1]+label
            curr_st = len(labels)*[int(curr_feat.location.start)-len(pf)] # get current start pos of motif
            curr_en = len(labels)*[int(curr_feat.location.end)-len(pf)] # get current end pos of motif
            annot_df = pd.concat([annot_df,pd.DataFrame([labels,curr_st,curr_en]).T]) # add info to initial df
    annot_df.columns = tf_loc_cols #change initial df column names
    annot_df.reset_index(drop=True, inplace=True) # reset index
    annot_df['TF'] = annot_df['TF'].str.upper() # change all TF names to uppercase
    lib_info = get_lib_info(lib) # get current lib info
    info_df = pd.DataFrame() #final df
    for i,pos in enumerate(lib_info['mut_loc']): # run over variable positions in library
        cond = (pos>=annot_df['Start']) & (pos<annot_df['End']) # find motifs in initial df containing the variable position
        sub_df = annot_df.loc[cond] # subset initial df by above condition
        for c_row in sub_df.index: # run over sub df
            row_info = list(sub_df.loc[c_row].values) #get info of row in initial df in list
            row_info.append(i) # add the current position
            info_df = pd.concat([info_df,pd.DataFrame(row_info).T]) # add TF motifs containing the variable position to the final df
    info_df.reset_index(drop=True, inplace=True) # reset index
    info_df.columns = tf_loc_cols + ['Pos'] # rename columns
    return info_df


if __name__ == '__main__':
    x = 1