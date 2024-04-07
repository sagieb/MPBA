import os
import numpy as np
import pandas as pd
import re
from functions import data_processing as dp
from functions import general_functions as gf
from functions import params
from collections import OrderedDict


def seq_mut_wt_cmp(res_table, lib_info, tf):
    '''This function gets a table, lib info and a tf name and returns the number of complete motifs
    of the tf found in each sequence of the table. in addition, it also returns the amount of variable
    positions that are wt in each sequence. the latter is not including the positions located in the tf motifs.
    the output format is a dataframe.'''
    tf_mot_loc = [i for i, s in enumerate(lib_info['mut_by_tf']) if re.split(r'(\d+)', tf)[0] in s] #  get poistions in
    # var nucs that are inside the tf motifs
    tf_mot_wt_seq = ''.join([lib_info['wt_at_var_loc'][i] for i in tf_mot_loc])  # get the wt seq within tf motifs
    # (all positions combined into string)
    wt_var_seq = ''.join(lib_info['wt_at_var_loc'])  # get wt seq and turn into a string
    wt_count = []
    wt_count_mot_pos = []
    for seq in res_table.index:  # iterating sequences from the input table
        wt_count.append(sum(map(str.__eq__, wt_var_seq, seq)))   # sum the number of wt nucleotides in
        # variable positions
        tf_mot_curr_seq = ''.join([seq[i] for i in tf_mot_loc])  # from the present sequence, take the nucs within tf
        # motif and make a string
        wt_count_mot_pos.append(sum(map(str.__eq__, tf_mot_wt_seq, tf_mot_curr_seq))) # sum the number of complete
        # tf motifs
    wt_count_norm = np.array(wt_count) - wt_count_mot_pos  # normalizing the wt count to the positions left
    # (ones without the tf motifs)
    wt_count_norm = list(wt_count_norm)
    oc_ic = pd.concat([pd.Series(wt_count_norm), pd.Series(wt_count_mot_pos)], axis=1)  # make the output dataframe
    oc_ic.index = res_table.index
    return oc_ic


def get_2_pos_seq_comb_id(pos_1_loc, pos_2_loc, wt_var_seq, res_table):
    '''This function get tow positions, for which it returns the 4 variables of
     indices of sequences in the library sequences:
    1. the first position is as in the wt variable sequence but the second position is mutated
    2. the second position is as in the wt variable sequence but the first position is mutated
    3. both positions are as in the wt variable sequence
    4.both positions are mutated
    the function gets as inputs the location of two positions within the variable sequence, the wt variable sequence and
     the library sequences.'''
    lib_seqs = list(res_table.index)
    pos_1_wt_nuc = wt_var_seq[pos_1_loc] # get the wt nucleotide in first position
    pos_2_wt_nuc = wt_var_seq[pos_2_loc] # get the wt nucleotide in second position
    pos1_wt_seqs_ids = [lib_seqs.index(seq) for seq in lib_seqs if (seq[pos_1_loc] == pos_1_wt_nuc)
                        & (seq[pos_2_loc] != pos_2_wt_nuc)]
    # get indices of sequences in which the first position is as in the wt variable but the second position is mutated
    pos2_wt_seqs_ids = [lib_seqs.index(seq) for seq in lib_seqs if (seq[pos_2_loc] == pos_2_wt_nuc)
                        & (seq[pos_1_loc] != pos_1_wt_nuc)]
    # get indices of sequences in which the second position is as in the wt variable sequence but
    # the first position is mutated
    both_wt_seqs_ids = [lib_seqs.index(seq) for seq in lib_seqs if (seq[pos_1_loc] == pos_1_wt_nuc)
                        & (seq[pos_2_loc] == pos_2_wt_nuc)]
    # get indices of sequences in which the both positions are as in the wt variable sequence
    both_mut_wt_seqs_ids = [lib_seqs.index(seq) for seq in lib_seqs if (seq[pos_1_loc] != pos_1_wt_nuc)
                            & (seq[pos_2_loc] != pos_2_wt_nuc)]
    # get indices of sequences in which the both positions are mutated
    return pos1_wt_seqs_ids, pos2_wt_seqs_ids, both_wt_seqs_ids, both_mut_wt_seqs_ids


def get_mot_mean_table(res_t, lib, tf):
    ''' this function takes a result table and generates an averaged table using only the positions localized in the
     motifs of the given tf (res_t_opt_mean) it also gives you the number of full tf motifs found in each averaged
     sequence (tf_mot_sum).'''
    curr_lib_info = dp.get_lib_info(lib)  # get library information
    wanted_ids = [i for i, s in enumerate(curr_lib_info['mut_by_tf']) if re.split(r'(\d+)', tf)[0] in s]  # get
    # positions of variable nucs found in the tf motifs
    sub_seqs = gf.list_substring(res_t.index, wanted_ids)  # get a list of substrings with only variable positions found
    # in tf motifs
    wt_subs_seq = gf.get_substring(curr_lib_info['wt_at_var_loc'], wanted_ids) # get wt nucleotides found in variable
    # positions within tf motifs
    full_mot_mat = np.zeros((len(sub_seqs), len(wanted_ids)))  # generate an empty matrix for binary intact motif in
    # each position of each sub sequence
    for i, seq in enumerate(sub_seqs):
        for pos in range(len(wanted_ids)):
            if seq[pos] == wt_subs_seq[pos]:
                   full_mot_mat[i,  pos] = 1  # 1 indicates that for this seqeuence in this position
                   # the tf motif was intact
    uniq_01, uniq_id = np.unique(full_mot_mat, axis=0, return_inverse=True)  # because we use only a subset of positions
    # a few sequence will contain the same sub sequence so we generate a unique matrix with all existing options
    # (uniq_01). in addition, we take a list telling us to which sub-sequence type each sequence belongs to (uniq_id)
    res_t_opt_mean = pd.DataFrame()
    for id in range(len(uniq_01)):  # running over all unique options
        curr_row = res_t.iloc[np.where(uniq_id == id)].mean().to_frame()  # averaging over all sequence with the same
        # sub-sequence of tf motif (intact/broken) combination and turning
        # into data frame
        curr_row = pd.pivot_table(curr_row, columns=curr_row.index, sort=False)  # using pivot table to switch the order
        # of the df to fit the following concatanation
        res_t_opt_mean = pd.concat([res_t_opt_mean, curr_row], axis=0)  # concatanating the new df into the existing df
        # containing all previous averaged sub-sequences
    id_names = []
    for patt in uniq_01: # generating names that will be used as index for the final table
        # (the format is 0_0_0/1_1_1 and all combinations in between)
        mid_list = []
        for charac in patt:
            mid_list.append(str(int(charac)))
        id_names.append('_'.join(mid_list))
    res_t_opt_mean.index = id_names  # changing the index names to the informative ones generated above
    tf_mot_sum = np.sum(uniq_01, axis=1)  # getting the number of full tf motifs found in each averaged sequence.
    return res_t_opt_mean, tf_mot_sum

def mean_vals_by_id(res_table, ids):
    '''This function returns the mean change (over bio repeats) relative to time point zero of specified sequences by
    tps The function gets the result table and the sequence ids.'''
    col_tps = [col.split('_') if len(col.split('_')) > 1 else col for col in res_table.columns.values]  #
    tps = list(OrderedDict.fromkeys(col_tps))  # get all time point by result table column names. convention: trans_tp
    pos_change = []  # empty list for holding tp average change
    for ti, tp in enumerate(tps):  # iterate over time points
        tp_col_ids = np.where(np.isin(col_tps, tp))[0]  # find the columns of the relevant time point
        pos_change.append(res_table.iloc[ids, tp_col_ids].mean().mean())  # calculate the mean
    return pos_change


def sem_vals_by_id(res_table, ids):
    '''This function returns the sem change (over bio repeats) relative to time point zero of specified sequences by
    tps The function gets the result table and the sequence ids.'''
    col_tps = [col.split('_') if len(col.split('_'))>1 else col for col in res_table.columns.values]  #
    tps = list(OrderedDict.fromkeys(col_tps)) # get all time point by result table column names. convention: trans_tp
    sem_over_ids = []  # empty list for holding tp average change
    for ti, tp in enumerate(tps):  # iterate over time points
        tp_col_ids = np.where(np.isin(col_tps, tp))[0]  # find the columns of the relevant time point
        sem_over_ids.append(res_table.iloc[ids, tp_col_ids].mean().sem())  # calculate the mean
    return sem_over_ids


def mean_pair_opts(pair, wt_var_seq, res_table):
    '''This function uses the functions get_2_pos_seq_comb_id and mean_vals_by_id and for a given pair of position
    and calculates the mean change relative to time point zero for three groups of sequences:
    1. the first position is as in the wt variable sequence but the second position is mutated
    2. the second position is as in the wt variable sequence but the first position is mutated
    3. both positions are as in the wt variable sequence
    the function gets as inputs a list of the location of two positions within the variable sequence, the wt variable
    sequence and the result table.'''
    pos_1 = pair[0]  # extract first position location
    pos_2 = pair[1]  # extract second position location
    pos1_wt_seqs_ids, pos2_wt_seqs_ids, both_wt_seqs_ids, _ = get_2_pos_seq_comb_id(pos_1, pos_2, wt_var_seq, res_table)
    # get the 3 group indices
    pos1_mean = mean_vals_by_id(res_table, pos1_wt_seqs_ids)
    pos2_mean = mean_vals_by_id(res_table, pos2_wt_seqs_ids)  # calculate the mean change for each group
    both_mean = mean_vals_by_id(res_table, both_wt_seqs_ids)
    return pos1_mean, pos2_mean, both_mean


def norm_non_cut(lib_num, res_table, avg_over=None, space=None):
    '''this function normalizes the result table to the mean of the number (avg_over) of sequences that were cut the
    least.Default parameters: avg_over=None - normalize to the mutant sequence Space: log2'''
    if avg_over==None:
            mut_var_seq = get_mut_var(lib_num)  # normalize to the mutant sequence
            if space == 'linear':
                norm_noncut = res_table / (res_table.loc[mut_var_seq])  # log2 space
            else:
                norm_noncut = res_table - (res_table.loc[mut_var_seq])  # linear space

    elif avg_over == 0:  # do not normalize to any sequence
        norm_noncut = res_table

    else: # normalize to n sequences that are the least cut
        if space == 'linear':  # this condition differentiate between linear and log space
            norm_noncut = res_table / \
                          (res_table.loc[res_table.mean(axis=1).sort_values(ascending=False)[0:avg_over].index]).mean()
        else:
            norm_noncut = res_table - \
                          (res_table.loc[res_table.mean(axis=1).sort_values(ascending=False)[0:avg_over].index]).mean()
    return norm_noncut


def get_only_pos_intact(relevant_locs, wt_var_seq, res_table):
    '''This function finds the sequence ids for which the given postions are intact and the rest are mutated.
    The function gets as inputs the a lists of position locations within the variable sequence, the wt variable
    sequence and the library sequences.'''
    rest_pos_locs = np.where(~np.isin(list(range(len(wt_var_seq))), relevant_locs))[0]  # find all position locations
    # except for the ones specified
    lib_seqs = list(res_table.index)
    wt_relevat_locs_nucs = gf.get_substring(wt_var_seq, relevant_locs)

    for i, seq in enumerate(lib_seqs):  # iterate over library sequences
        broken_pos_num = 0  # initialize a variable for counting the the number of mutated positions
        seq_relevat_locs_nucs = gf.get_substring(seq, relevant_locs)
        only_pos_intact_ids = []
        for pos in rest_pos_locs:  # iterate over positions locations of the not specified positions
            if (seq_relevat_locs_nucs == wt_relevat_locs_nucs) & (seq[pos] != wt_var_seq[pos]):
                broken_pos_num += 1
                # if the sequnce contains the wt option in the specified position but a mutation in the other position
                # update broken_pos_num
        if len(rest_pos_locs) == broken_pos_num:
            # in case the all of the positions are mutated except for the specified position, return the sequence index
            only_pos_intact_ids = i
            break
    return only_pos_intact_ids


def get_only_pos_mut(relevant_seqs, lib_num):
    '''This function gets a list of the res_table relevant sequences and the library number and
    returns a list of indices of the sequences in which only one position is mutated and the rest are intact.
    Mutation positions correspond the index position in the output list.'''
    lib_info_dict = dp.get_lib_info(lib_num)  # get lib info
    wt_var_seq = ''.join(lib_info_dict['wt_at_var_loc'])  # get wt sequence
    pos_opts = lib_info_dict['var_pos_pattern']  # get a string of all options per position (for example:'[AG][AG]')
    only_pos_mut_seq_ids = []  # initialize a variable for storing the indices for sequences with only one mutation
    for i, pos in enumerate(wt_var_seq):  # iterate the wild type sequence positions
        curr_opst = np.array(list(re.sub(r'[^a-zA-Z]', '', pos_opts.split('][')[i]).strip()))
        # get an array of current position options
        mut_opt = curr_opst[curr_opst != pos]  # get the mutation option
        wt_var_seq_list = list(wt_var_seq)  # replace the wt position in the mutated position
        wt_var_seq_list[i] = mut_opt[0]
        only_pos_mut_seq = ''.join(wt_var_seq_list)
        only_pos_mut_seq_ids.append(relevant_seqs.index(only_pos_mut_seq))
        # find the mutated sequence in the relevant sequences list and store in only_pos_mut_seq_ids
    return only_pos_mut_seq_ids

def convert_fc_to_occ(vals):
    '''This function gets log2 fold change values / dataframe and returns the calculated occupancy.'''
    frac_cut = 1 - (0.5 ** (-vals))
    return frac_cut

def norm_res_data(tf, lib_num, exp_num, norm_to, occ=None):
    '''This function gets the tf name, library number and experiment number,
    and returns a normalized data frame by these steps:
    1. norm_reads  2. remove bad samples  3. mean_bio_rep_tps  4.res_log2  5.norm_to_tp_0
    6. norm_non_cut 7.convert to occupancy (optional).
    optional parameter occ: if true: returns the percentage occupancy, default returns log2 fold change.'''
    file_name = tf + '_' + str(lib_num) + '_' + str(exp_num) + '.csv'  # find results file ('tf_libnum_exp.csv')
    res_table = (pd.read_csv(os.path.join(params.RES_PATH, file_name), index_col=0))
    norm_res = dp.norm_reads(res_table)  # 1
    sample_filt_norm, _ = dp.rm_samples(tf, lib_num, exp_num, norm_res)  # 2
    sample_info = dp.get_samp_info(sample_filt_norm)
    mean_repeats = dp.mean_bio_rep_tps(sample_filt_norm, sample_info)  # 3
    log2_norm = dp.res_log2(mean_repeats)  # 4
    norm_tp_0 = dp.norm_to_tp_0(log2_norm, sample_info)  # 5
    norm_uncut = norm_non_cut(lib_num, norm_tp_0, norm_to)  # 6
    if occ == True:  # 7
        norm_occ = convert_fc_to_occ(norm_uncut)
        return norm_occ
    else:
        return norm_uncut


def norm_without_biorep_mean(tf, lib_num, exp_num, norm_to, occ=None):
    file_name = tf + '_' + str(lib_num) + '_' + str(exp_num) + '.csv'  # find results file ('tf_libnum_exp.csv')
    res_table = (pd.read_csv(os.path.join(params.RES_PATH, file_name), index_col=0))
    norm_res = dp.norm_reads(res_table)  # 1
    sample_filt_norm, _ = dp.rm_samples(tf, lib_num, exp_num, norm_res)  
    sample_info = dp.get_samp_info(sample_filt_norm)
    log2_norm = dp.res_log2(sample_filt_norm)  
    norm_tp_0 = dp.norm_to_tp_0(log2_norm, sample_info)  
    norm_uncut = norm_non_cut(lib_num, norm_tp_0, norm_to) 
    if occ == True:  # 7
        norm_occ = convert_fc_to_occ(norm_uncut)
        return norm_occ
    else:
        return norm_uncut


def mean_pos_change_calc(tf, lib_num, exp_num, norm_to, occ=None):
    '''This function calculates the mean position change, by averaging over all sequences in which a position is intact.
    The function return a data frame containing the mean position change (tps X positions).
    when occ is true the returned values will by in occupancy.'''
    norm_df = norm_res_data(tf, lib_num, exp_num, norm_to, occ)  # normalize data
    lib_info_dict = dp.get_lib_info(lib_num)  # get library info
    wt_var_seq = ''.join(lib_info_dict['wt_at_var_loc'])  # get wild-type sequence
    relevant_seqs = list(norm_df.index)  # get a list of all variants
    col_tps = [x.split('_')[1] for x in norm_df.columns.values]  # get experiment tps
    tps = list(OrderedDict.fromkeys(col_tps))
    mean_pos_change = pd.DataFrame(index=tps)  # initialize a variable for storing the mean position change
    # (tps X positions)
    for i,curr_pos in enumerate(wt_var_seq):  # iterate over wild-type sequence positions
        relevant_seqs_ids = [j for j,seq in enumerate(relevant_seqs) if seq[i] == curr_pos]
        # get the sequence indices in which the position i is wt
        curr_pos_change = mean_vals_by_id(norm_df, relevant_seqs_ids)
        # mean over the wt-position sequence indices found
        mean_pos_change[str(i)] = curr_pos_change  # add the position mean change to the mean_pos_change data frame
    mean_pos_change.columns = lib_info_dict['mut_by_tf']  # set column named as the TF binding each position
    return mean_pos_change


def one_pos_mut_rest_mean_clac(tf, lib_num, exp_num, norm_to, occ=None):
    '''his function calculates the mean position change, by averaging over all sequences in which a position is mutated
    The function return a data frame containing the mean position change (tps X positions).
    when occ is true the returned values will by in occupancy.'''
    norm_df = norm_res_data(tf, lib_num, exp_num, norm_to, occ)  # normalize data
    lib_info_dict = dp.get_lib_info(lib_num)  # get library info
    wt_var_seq = ''.join(lib_info_dict['wt_at_var_loc'])  # get wild-type sequence
    col_tps = [x.split('_')[1] for x in norm_df.columns.values]  # get experiment tps
    tps = list(OrderedDict.fromkeys(col_tps))
    seqs = list(norm_df.index)  # get a list of all variants
    pos_opts = lib_info_dict['var_pos_pattern']
    mut_var = ''  # Generate the sequences that is fully mutated
    for i, pos in enumerate(wt_var_seq):  # iterate the wild type sequence positions
        curr_opst = np.array(list(re.sub(r'[^a-zA-Z]', '', pos_opts.split('][')[i]).strip()))
        mut_var += curr_opst[curr_opst != wt_var_seq[i]][0]

    mut_pos_mean = pd.DataFrame(index=tps, columns=[])
    for i, pos in enumerate(mut_var):
        seq_ids = [si for si, s in enumerate(seqs) if s[i] == pos]
        mut_pos_mean[str(i)] = mean_vals_by_id(norm_df, seq_ids)
    mut_pos_mean.columns = lib_info_dict['mut_by_tf']  # set column names as the TF binding each position
    return mut_pos_mean


def only_pos_intact_change_calc(tf, lib_num, exp_num=None, norm_to=None, occ=None, norm_df=None):
    '''This function calculates for each position the change relative to time point 0 for the single sequence in which
    the postion is wild-type and the rest are mutates. The function returns a data frame containing these changes for
    each position (tps X positions).
    IMPORTANT: this function returns the mean over biological repeats!.'''
    if norm_df is None:
        norm_df = norm_res_data(tf, lib_num, exp_num, norm_to, occ)  # normalize data
    lib_info_dict = dp.get_lib_info(lib_num)  # get library info
    wt_var_seq = ''.join(lib_info_dict['wt_at_var_loc'])  # get wild-type sequence
    col_tps = [x.split('_')[1] for x in norm_df.columns.values]  # get experiment tps
    tps = list(OrderedDict.fromkeys(col_tps))
    only_intact_pos_change = pd.DataFrame(index=tps)
    only_intact_pos_change_sem = pd.DataFrame(index=tps)
    # initialize a variable for storing the sequence change, when only one position is wt (tps X positions)
    for i in range(len(wt_var_seq)):  # iterate over wild-type sequence positions
        only_pos_seq_id = get_only_pos_intact([i], wt_var_seq, norm_df)  # get the sequence id for which only the
        # position i is wild-type
        mean_seq_tps = dp.mean_over_bio_reps(norm_df.iloc[only_pos_seq_id, :].to_frame().transpose()).values[0]
        # get the fold change values of the sequence found
        only_intact_pos_change[str(i)] = mean_seq_tps
        # add the sequence change to the only_intact_pos_change data frame

        sem_seq_tps = dp.sem_over_bio_reps(norm_df.iloc[only_pos_seq_id, :].to_frame().transpose()).values[0]
        # get the fold change values of the sequence found
        only_intact_pos_change_sem[str(i)] = sem_seq_tps  # add the sequence change to the only_intact_pos_change df

    only_intact_pos_change.columns = lib_info_dict['mut_by_tf']  # set column names as the TF binding each position
    only_intact_pos_change_sem.columns = lib_info_dict['mut_by_tf']
    return only_intact_pos_change, only_intact_pos_change_sem


def only_pos_mut_change_calc(tf, lib_num, exp_num, norm_to, occ=None):
    '''This function calculates for each position the change relative to time point 0 for the single sequence in which
    the position is mutated and the rest are wild-type. The function returns a data frame containing these changes for
    each position (tps X positions).'''
    norm_df = norm_res_data(tf, lib_num, exp_num, norm_to, occ)  # normalize data
    lib_info_dict = dp.get_lib_info(lib_num)  # get library info
    col_tps = [col.split('_')[1] for col in norm_df.columns.values]  # get experiment tps
    tps = list(OrderedDict.fromkeys(col_tps))
    relevant_seqs = list(norm_df.index)  # get a list of all variants
    one_mut_ids = get_only_pos_mut(relevant_seqs, lib_num)
    # get a list of indices of the sequences in which only one position is mutated and the rest are intact
    only_mut_pos_change = pd.DataFrame(index=tps)
    # initialize a variable for storing the sequence change, when only one position is mutated (tps X positions)
    for i, curr_id in enumerate(one_mut_ids):  # iterate over the ids of sequences containing a mutation
        mean_seq_tps = dp.mean_over_bio_reps(norm_df.iloc[curr_id, :].to_frame().transpose()).values[0]
        # get the fold change values of the sequence found
        only_mut_pos_change[str(i)] = mean_seq_tps  # add the sequence change to the only_mut_pos_change data frame
    only_mut_pos_change.columns = lib_info_dict['mut_by_tf']  # set column names as the TF binding each position
    return only_mut_pos_change


def generate_pos_change_df(tf, lib_num, exp_num, norm_to, occ=None):
    '''This function returns a data frame that contains:
    the mean values when each position is intact and the other is not, when both are wt, when both mutated.
    in addition the dataframe contains the single sequence in which both positions are wt and both positions
    are mutated.'''
    norm_df = norm_res_data(tf, lib_num, exp_num, norm_to, occ)  # normalize data
    lib_info_dict = dp.get_lib_info(lib_num)  # get library info
    wt_var_seq = ''.join(lib_info_dict['wt_at_var_loc'])  # get wild-type sequence
    col_tps = [col.split('_')[1] for col in norm_df.columns.values]  # get experiment tps
    tps = list(OrderedDict.fromkeys(col_tps))
    pos_num = len(wt_var_seq)  # number of variable positions
    pair_combs = gf.get_pair_combs(list(range(pos_num)), 2)  # get all pairwise combination of variable positions
    seqs = list(norm_df.index)

    col_names = ['Pos1', 'Pos2', 'Time_point', 'Pos1_only_mean', 'Pos2_only_mean', 'Both_wt_mean', 'Both_mut_mean',
                 'Both_wt', 'Both_mut']
    comb_df = pd.DataFrame(columns=col_names)
    # initialize a variable for storing fold change (relative to tp0) by position pair in a certain time-point
    for comb in pair_combs:
        pos_1_loc = comb[0]  # position 1
        pos_2_loc = comb[1]  # position 2
        only_pos1_wt_seqs_ids, only_pos2_wt_seqs_ids, both_wt_seqs_ids, both_mut_seqs_ids = get_2_pos_seq_comb_id(
            pos_1_loc, pos_2_loc, wt_var_seq, norm_df)
        # get sequences id groups for the specified position: 1. pos1 wt pos2 mut, 2. pos1 wt pos2 xwt,
        # 3. both wt, 4. both mut
        pos_1_only_mean_vals = mean_vals_by_id(norm_df, only_pos1_wt_seqs_ids)
        pos_2_only_mean_vals = mean_vals_by_id(norm_df, only_pos2_wt_seqs_ids)
        both_wt_mean_vals = mean_vals_by_id(norm_df, both_wt_seqs_ids)  # mean the change by id groups
        both_mut_mean_vals = mean_vals_by_id(norm_df, both_mut_seqs_ids)
        pos_opts = lib_info_dict['var_pos_pattern']
        mut_var = ''  # Generate the sequences that is fully mutated
        for i, pos in enumerate(wt_var_seq):  # iterate the wild type sequence positions
            curr_opst = np.array(list(re.sub(r'[^a-zA-Z]', '', pos_opts.split('][')[i]).strip()))
            mut_var += curr_opst[curr_opst != wt_var_seq[i]][0]

        pair_intact_only_seq = list(mut_var)
        # generate the sequence in which only the pair is wt and the rest are mutated
        pair_intact_only_seq[pos_1_loc] = wt_var_seq[pos_1_loc]
        pair_intact_only_seq[pos_2_loc] = wt_var_seq[pos_2_loc]
        pair_intact_only_seq = "".join(pair_intact_only_seq)

        pair_mut_only_seq = list(wt_var_seq)
        pair_mut_only_seq[pos_1_loc] = mut_var[pos_1_loc]
        pair_mut_only_seq[pos_2_loc] = mut_var[pos_2_loc]
        pair_mut_only_seq = "".join(pair_mut_only_seq)

        only_comb_wt_seq_id = seqs.index(pair_intact_only_seq)  # get the fold change values by the ids found
        only_comb_mut_seq_id = seqs.index(pair_mut_only_seq)
        only_comb_wt_vals = mean_vals_by_id(norm_df, [only_comb_wt_seq_id])
        only_comb_mut_vals = mean_vals_by_id(norm_df, [only_comb_mut_seq_id])

        for i, tp in enumerate(tps):  # iterate over time points
            comb_tp_vals = pd.DataFrame([[pos_1_loc, pos_2_loc, tp, pos_1_only_mean_vals[i], pos_2_only_mean_vals[i],
                                          both_wt_mean_vals[i], both_mut_mean_vals[i], only_comb_wt_vals[i],
                                          only_comb_mut_vals[i]]], columns=col_names)
            comb_df = pd.concat([comb_df, comb_tp_vals])
            # update comb_df

    for tp in tps:  # iterate over time points
        for pi in range(pos_num):  # add self combination to the data frame with 0 values
            comb_df = pd.concat([comb_df, pd.DataFrame([[pi, pi, tp, 0, 0, 0, 0, 0, 0]], columns=col_names)])
    return comb_df


def get_extreme_cut(tf, lib_num, exp, num_seqs):
    '''This function finds the indexes of the most and least cut sequences.
    The number of sequences you'll get is defined by `num_seqs`
    It has two outputs - the first are the indexes of the most cut sequences
    and the second is of the least cut'''
    file_name = tf + '_' + str(lib_num) + '_' + str(exp) + '.csv'  # find results file ('tf_libnum_exp.csv')
    res_table = (pd.read_csv(os.path.join(params.RES_PATH, file_name), index_col=0))
    norm_data = dp.norm_reads(res_table)  # 1
    filt_norm, _ = dp.rm_samples(tf, lib_num, exp, norm_data)
    sample_info = dp.get_samp_info(filt_norm)
    norm_mean_tech = dp.mean_bio_rep_tps(filt_norm, sample_info).apply(np.log2)
    sample_info = dp.get_samp_info(norm_mean_tech)
    relative_tp0_norm = dp.norm_to_tp_0(norm_mean_tech, sample_info)
    most_cut_id = relative_tp0_norm.median(axis=1).sort_values().index[0:num_seqs]
    least_cut_id = relative_tp0_norm.median(axis=1).sort_values(ascending=False).index[0:num_seqs]
    return most_cut_id, least_cut_id


def get_all_tf_motifs_wt_or_mut(tf, lib_info, lib_seqs):
    lib_seqs = list(lib_seqs)
    wt_var_seq = ''.join(lib_info['wt_at_var_loc'])  # get wild-type sequence
    tf_name_split = re.split(r'(\d+)', tf)[0]
    tf_pos = [i for i, mot in enumerate(lib_info['mut_by_tf']) if tf_name_split in mot]
    non_tf_pos = np.setdiff1d(range(len(wt_var_seq)), tf_pos)

    mut_pos = {}
    intact_pos = {}
    for pos in range(len(wt_var_seq)):
        intact_pos[pos] = [lib_seqs.index(seq) for seq in lib_seqs if (seq[pos]==wt_var_seq[pos])]
        mut_pos[pos] = [lib_seqs.index(seq) for seq in lib_seqs if (seq[pos]!=wt_var_seq[pos])]

    tf_pos_wt_ids = list(range(len(lib_seqs)))
    tf_pos_mut_ids = list(range(len(lib_seqs)))
    for pos in tf_pos:
        tf_pos_wt_ids = np.intersect1d(tf_pos_wt_ids, intact_pos[pos])
        tf_pos_mut_ids = np.intersect1d(tf_pos_mut_ids, mut_pos[pos])

    non_tf_pos_wt_ids = list(range(len(lib_seqs)))
    non_tf_pos_mut_ids = list(range(len(lib_seqs)))
    for pos in non_tf_pos:
        non_tf_pos_wt_ids = np.intersect1d(non_tf_pos_wt_ids, intact_pos[pos])
        non_tf_pos_mut_ids = np.intersect1d(non_tf_pos_mut_ids, mut_pos[pos])

    return tf_pos_wt_ids, tf_pos_mut_ids, intact_pos, mut_pos, non_tf_pos_wt_ids, non_tf_pos_mut_ids


def get_mut_var(lib_num):
    lib_info_dict = dp.get_lib_info(lib_num)  # get lib info
    wt_var_seq = ''.join(lib_info_dict['wt_at_var_loc'])  # get wt sequence
    pos_opts = lib_info_dict['var_pos_pattern']  # get a string of all options per position (for example:'[AG][AG]')

    mut_seq = ''  # initialize a variable for storing the indices for sequences with only one mutation
    for i, pos in enumerate(wt_var_seq):  # iterate the wild type sequence positions
        curr_opst = np.array(list(re.sub(r'[^a-zA-Z]', '', pos_opts.split('][')[i]).strip()))  # get an array of
        # current position options
        mut_seq = mut_seq + curr_opst[curr_opst != pos][0]  # get the mutation option
    return mut_seq


def get_tf_mot_pos(tf, lib_num):
    '''This function returns a list containing the positions of the tf motifs in the wanted library'''
    lib_info = dp.get_lib_info(lib_num)  # get lib info
    tf_red = re.findall('\d*\D+',tf)[0]  # take tf letters without number to fit lib info
    # (for duplicated we wrote Nrg and not Nrg1/Nrg2)
    split_by_pos = [pos.split('/') for pos in lib_info['mut_by_tf']]  # list of lists. each sublist contain all tf
    # motifs mutated in this position
    tf_mot_pos = [i for i,curr_pos in enumerate(split_by_pos) if
                  any(tf_red in curr_tf_in_pos for curr_tf_in_pos in curr_pos)] # get position of tf motifs
    return tf_mot_pos


def non_tf_mot_pos(tf, lib_num):
    '''This function returns all non tf motif positions in the library of interest.'''
    lib_info = dp.get_lib_info(lib_num)  # get lib info
    all_pos_list = np.arange(0, len(lib_info['wt_at_var_loc']))  # an array of all positions of the library
    tf_pos_list = get_tf_mot_pos(tf, lib_num)  # get tf motifs positions
    nontf_pos_list = np.setdiff1d(all_pos_list, tf_pos_list)  # find non tf motif positions
    return nontf_pos_list


def get_values_by_wt_pos(pos_list, norm_df, lib_num, mean_seqs=None):
    '''This function gets position list, return the value of the case in which the specifies positions are wild-type.
    default: all the rest of the positions are mutated.
    when mean_seqs=True: mean over all other positions.'''
    seqs = list(norm_df.index)
    all_seqs_split = np.array([list(seq) for seq in seqs])
    lib_info = dp.get_lib_info(lib_num)  # get library info
    wt_var_seq = ''.join(lib_info['wt_at_var_loc'])  # get wild-type sequence
    mut_var_seq = get_mut_var(lib_num)  # get mutates sequence
    non_pos = np.delete(np.array(range(0, len(wt_var_seq))), pos_list)  # all positions other than relevant positions
    wt_seq_in_pos = np.array(list(wt_var_seq))[pos_list]  # wt seq at relevant positions
    mut_seq_in_non_pos = np.array(list(mut_var_seq))[non_pos]  # mut seq at non-relevant positions
    id_rel_wt = np.sum(all_seqs_split[:, pos_list] == wt_seq_in_pos, axis=1) == len(pos_list)
    id_rel_mut = np.sum(all_seqs_split[:, non_pos] == mut_seq_in_non_pos, axis=1) == len(non_pos)
    # ids of seqs with mut seq at non-relevant positions
    if mean_seqs==None:
        ids = id_rel_wt & id_rel_mut
        seq_values = norm_df.loc[np.array(seqs)[ids], :]
    else:
        ids = id_rel_wt
        seq_values = norm_df.loc[seqs[id_rel_wt]].mean()
    ids = np.where(ids)[0]
    return seq_values, ids


def get_values_by_mut_pos(pos_list, norm_df, lib_num, mean_seqs=None):
    '''This function gets position list, return the value of the case in which the specifies positions are mutated.
    default: all the rest of the positions are wild-type.
    when mean_seqs=True: mean over all other positions.'''
    seqs = list(norm_df.index)
    all_seqs_split = np.array([list(seq) for seq in seqs])
    lib_info = dp.get_lib_info(lib_num)  # get library info
    wt_var_seq = ''.join(lib_info['wt_at_var_loc'])  # get wild-type sequence
    mut_var_seq = get_mut_var(lib_num)  # get mutates sequence
    non_pos = np.delete(np.array(range(0, len(wt_var_seq))), pos_list)  # all positions other than relevant positions
    wt_seq_in_pos = np.array(list(wt_var_seq))[non_pos]  # wt seq at all positions except for the specified one
    mut_seq_in_non_pos = np.array(list(mut_var_seq))[pos_list]  # mut seq at specified positions
    id_wt_pos = np.sum(all_seqs_split[:, non_pos] == wt_seq_in_pos, axis=1) == len(non_pos)
    id_mut_pos = np.sum(all_seqs_split[:, pos_list] == mut_seq_in_non_pos, axis=1) == len(pos_list)
    # ids of seqs with mut seq in specified positions
    if mean_seqs == None:
        seq_values = norm_df.loc[np.array(seqs)[id_wt_pos & id_mut_pos], :]
    else:
        seq_values = norm_df.loc[np.array(seqs)[id_mut_pos]].mean()

    return seq_values


def obs_calc(tf, lib_num, exp_num, relevant_pos, relevant_combs, tp, avg=False, direction='gain'):
    '''This function calculates the observed values in three methods:
    (1) "single gain" - the contribution of the relevant positions when observing a single sequence in which only
        the specified positions are wt.
    (2) "mean gain" - the contribution of the relevant positions when calculating the mean over all sequences in which
        the specified positions are wt, normalized to the mean of sequences in which the specified positions are
        mutated.
    (3) "single loss - the effect from the loss of the relevant positions when observing the sequence in which only the
        position is mutated. This is the difference between the wild-type sequence value and the sequence in which the
        specified positions are mutated and the rest are wild-type.
    The output is a dictionary in which the keys are the positions and values are the observed calculations
    in occupancy.'''
    lib_info = dp.get_lib_info(lib_num)  # get library info
    wt_var_seq = ''.join(lib_info['wt_at_var_loc'])  # get wild-type sequence
    observed_dict = {}  # initialize an empty dictionary to store the observed values by positions (keys)

    if avg == False and direction == 'gain':  # calculate according to first method
        norm_df = dp.mean_over_bio_reps(
        norm_res_data(tf, lib_num, exp_num, norm_to=None))
        # get the mean normalized data in fold change values, normalize to mutant
        only_intact_pos_change = only_pos_intact_change_calc(tf, lib_num, exp_num, norm_to=None)  # get dataframe of
        # positionsX time points, values are the log2 fold change of sequences in which only the position is wt.
        only_intact_pos_change = only_intact_pos_change.iloc[only_intact_pos_change.index == tp, :]
        # only_intact_pos_change by specified time point
        for pos in relevant_pos:  # iterate over positions
            observed_dict[str(pos)] = convert_fc_to_occ(only_intact_pos_change.values[0][pos])
            # get values for single positions
        for comb in relevant_combs: # iterate over combinations
            seq_val = get_values_by_wt_pos(comb, norm_df, lib_num)[tp][0]
            # get the value for the sequence in which only the combination is wt
            observed_dict[''.join(map(str, comb))] = convert_fc_to_occ(seq_val)

    elif avg == True and direction == 'gain': # calculate according to second method
        norm_df = dp.mean_over_bio_reps(
        norm_res_data(tf, lib_num, exp_num, norm_to=0))  # normalize data without the normalization to mutant seq
        mean_comb_df = generate_pos_change_df(tf, lib_num, exp_num, norm_to=0)  # get combination calculations
        for comb in relevant_combs:  # iterate over combinations
            pos1 = comb[0]  # first position
            pos2 = comb[1]  # second position
            # pos_list = np.delete(np.array(range(0, len(wt_var_seq))), comb)
            # get all position indices except for the combination
            # comb_mut_val = lto.get_values_by_wt_pos(pos_list, norm_df, lib_num)[tp][0]
            #########################################################################
            comb_mut_val = get_values_by_mut_pos(comb, norm_df, lib_num, mean_seqs=True)[tp]
            # get the value of the sequence in which the combination is mutated and mean overt the other positions are
            curr_comb_df = mean_comb_df.loc[(mean_comb_df['Pos1'] == pos1) &
                                            (mean_comb_df['Pos2'] == pos2) &
                                            (mean_comb_df['Time_point'] == tp)]
            # get the relevant row in the combination df
            observed_dict[str(pos1)] = convert_fc_to_occ(curr_comb_df['Pos1_only_mean'][0] - comb_mut_val)
            # (pos1 wt, pos2 mut, other mean) - (pos1+2 mut, other mean )
            observed_dict[str(pos2)] = convert_fc_to_occ(curr_comb_df['Pos2_only_mean'][0] - comb_mut_val)
            # (pos2 wt, pos1 mut, other mean) - (pos1+2 mut, other mean )
            observed_dict[str(pos1)+str(pos2)] = convert_fc_to_occ(curr_comb_df['Both_wt_mean'][0] - comb_mut_val)
            # (pos1+2 wt) - (pos1+2 mut, other mean )

    elif avg == False and direction == 'loss':
        norm_df = dp.mean_over_bio_reps(norm_res_data(tf, lib_num, exp_num, norm_to=0))
        # normalize data without the normzlization to mutant seq
        only_mut_pos_change = only_pos_mut_change_calc(tf, lib_num, exp_num, norm_to=0)
        # get dataframe of positionsX time points, values are the log2 fold change of sequences in which only the
        # position is mut and rest are wt.
        only_mut_pos_change = only_mut_pos_change.loc[tp].to_frame().T  # only_mut_pos_change in specified time point
        wt_var_seq_val = norm_df.loc[wt_var_seq][tp]  # wt sequence value
        for pos in relevant_pos:  # iterate over positions
            observed_dict[str(pos)] = convert_fc_to_occ(wt_var_seq_val - only_mut_pos_change.values[0][pos])
            # (wt seq) - (pos mut other wt)
        for comb in relevant_combs:  # iterate over combinations
            pos_list = np.delete(np.array(range(0, len(wt_var_seq))), comb)  # position indices of non combination
            seq_val = get_values_by_wt_pos(pos_list, norm_df, lib_num)[tp][0]  # comb mut, rest wt
            observed_dict[''.join(map(str, comb))] = convert_fc_to_occ(wt_var_seq_val - seq_val)
            # (wt seq) - (comb mut, rest wt)
    return observed_dict


def exp_calc(relevant_combs, observed_dict):
    '''This function uses the output dictionary of the function "obs_calc" and calculates the expected values for
    independency. Values are in occupancy.'''
    expected_dict = {}
    for comb in relevant_combs:   # iterate over combinations
        sub_comb = gf.get_pair_combs(comb, len(comb) - 1) # get all the sub combinations of current combination
        sub_comb_observed = []
        for curr_sub_comb in sub_comb:  # iterate over sub-combinations
            curr_key = ''.join(map(str, curr_sub_comb))  # generate key name
            sub_comb_observed.append(observed_dict[curr_key])  # get sub-combination observed value, and store in
            # sub_comb_observed list
        expected_dict[''.join(map(str, comb))] = (1 - np.prod(1 - (np.array(sub_comb_observed))))
        # expected calculation
    return expected_dict


def obs_exp_mean_context(tf, lib_num, tp, relevant_tfs, norm_df=None, exp_num=None):
    col_names_mean = ['Gene', 'Lib_num', 'Comb', 'Min_dist', 'Mut_dist', 'Comb_tfs', 'Observed_mean',
                      'Expected_mean', 'Obs-Exp_mean', 'Obs-Exp_sem', 'TF']
    mean_comb_df = pd.DataFrame(columns=col_names_mean)
    all_combs_all_vals = pd.DataFrame()

    lib_info = dp.get_lib_info(lib_num) # get library information
    wt_var_seq = ''.join(lib_info['wt_at_var_loc'])  # wt seq
    mut_var_seq = get_mut_var(lib_num) # mutated seq
    pos_num = len(wt_var_seq)  # number of variable positions in library
    all_pos = list(range(pos_num))  # all positions ids
    comb_size = 2 # combination size
    all_combs = gf.get_pair_combs(all_pos, 2)  # all possible combinations i library
    all_combs = np.reshape(all_combs, (len(all_combs), comb_size))  # all combination matrix (rows combinations)

    split_tfs = [re.match(r"([a-z]+)([0-9]+)", curr_tf, re.I).groups()[0] for curr_tf in relevant_tfs]
    # TF names without numbers
    tf_pos = []
    for split_tf in split_tfs: # iterate over tfs
        tf_pos += [i for i, pos in enumerate(lib_info['mut_by_tf']) if split_tf in pos]
        # store tf variable positions in "tf_pos"

    tf_combs_ids = np.unique(np.concatenate([np.where(np.sum(all_combs == pos, axis=1))[0] for pos in tf_pos]))
    #iterate over "tf_pos" and get all the ids of combinations containing this positions
    relevant_combs = all_combs[tf_combs_ids, :]  # get the combinations containing tf positions

    if norm_df is None:  # if the normalized result table was not provided, get it by input library information
        norm_df = dp.mean_over_bio_reps(norm_res_data(tf, lib_num, exp_num, norm_to=0))  # get the mean normalized
        # data in fold change values
    seqs = norm_df.index  # library sequences
    all_seqs_split = np.array([list(seq) for seq in seqs])  # library sequences split

    for comb in relevant_combs:  # iterate over combinations
        pos1 = comb[0]  # combination 1
        pos2 = comb[1]  # combination 2
        pos_dist = abs(lib_info['mut_loc'][pos1] - lib_info['mut_loc'][pos2])  # mutation distance
        comb_tfs = lib_info['mut_by_tf'][pos1] + '_' + lib_info['mut_by_tf'][pos2]  # TFs in combination
        split_comb_tfs = comb_tfs.split('_')  # tfs relevant to the combination

        tf_pos_i = [pos_i for pos_i, tfs in enumerate(split_comb_tfs) if tf[0:3] in tfs]
        # find the relevant tf position
        if len(tf_pos_i) == 0:
            min_dist = None
        else:
            tf_pos_i = tf_pos_i[0]
            second_tf_pos_i = abs(tf_pos_i - 1)
            pos_by_tf = np.array(['xxxx'] * len(comb))
            pos_by_tf[tf_pos_i] = tf
            if '/' in split_comb_tfs[second_tf_pos_i]:
                second_mot_tfs = split_comb_tfs[second_tf_pos_i].split('/')
                min_dist = 1000
                for curr_tf in second_mot_tfs:
                    pos_by_tf[second_tf_pos_i] = curr_tf
                    curr_dist = mot_dist_by_mut_comb(lib_num, comb, pos_by_tf)
                    min_dist = min(min_dist, curr_dist)
            else:
                pos_by_tf[second_tf_pos_i] = split_comb_tfs[second_tf_pos_i]
                min_dist = mot_dist_by_mut_comb(lib_num, comb, pos_by_tf)

        comb_name = str(pos1) + '_' + str(pos2)  # combination name
        comb_rm_seqs = np.delete(all_seqs_split, np.array(comb), axis=1)  # remove combination from "all_seqs_split"
        flanking_pos_seq_opts = np.unique([''.join(seq) for seq in comb_rm_seqs])  # join all remaining non-comb
        # positions to string sequences
        wt_opt = list(np.array(list(wt_var_seq))[np.array(comb)])  # get the wt option for combination
        mut_opt = list(np.array(list(mut_var_seq))[np.array(comb)])  # get the mutated option for combination
        comb_pos_seq = [[mut_var_seq[pos], wt_var_seq[pos]] for pos in comb]
        # generate a list containing for each position a sublist os the mutated option and the wt option
        comb_opts = gf.comb_2_arrays(comb_pos_seq[0], comb_pos_seq[1])  # generate all permutations of position
        # option sequence (++, --, +-, -+)
        mut_opt_i = [i for i, opt in enumerate(comb_opts) if list(opt) == mut_opt][0]
        # get the id of the mutated option in "comb_opta"
        comb_opts.pop(mut_opt_i) # remove the mutation option from "comb_opts"

        # order opts
        wt_opt_i = [i for i, opt in enumerate(comb_opts) if list(opt) == wt_opt][0]
        # get the id of the wt option in "comb_opta"
        comb_opts_ids = list(range(len(comb_opts)))  # the number of possible combination without the mutated sequence
        comb_opts_ids.pop(wt_opt_i)  # remove the wt option id from "comb_opts_ids"
        opts_order = [wt_opt_i, 0, 0] # fix position order (++, +-, -+)
        curr_opt = comb_opts[comb_opts_ids[0]]  # get the first combination opts that are not wt or mutated
        if curr_opt[0] == wt_opt[0]:  # if the first position in the option is wt
            opts_order[1] = comb_opts_ids[0]  # this option will be the second item in opts_order
            opts_order[2] = comb_opts_ids[1]  # and the other option will be the third
        else:  # if the first position in the option is mutated
            opts_order[2] = comb_opts_ids[0]  # this option will be the third item in opts_order
            opts_order[1] = comb_opts_ids[1]  # and the other option will be the second
        comb_opts = np.array(comb_opts)[np.array(opts_order)]  # order comb_opts by opts_order

        comb_df = pd.DataFrame()  # generate an empty data frame for the current combination
        for opt in comb_opts:  # iterate over options
            opt_vals = []
            contexts = []
            for seq in flanking_pos_seq_opts:  # iterate over context position sequences
                curr_seq = list(seq)
                mut_seq = list(seq)
                for pos in range(comb_size):  # iterate over position ids in combination
                    curr_seq.insert(comb[pos], opt[pos])  # insert current option in current context
                    mut_seq.insert(comb[pos], mut_opt[pos])  # insert mutated option in current context
                curr_seq = ''.join(curr_seq)
                mut_seq = ''.join(mut_seq)
                opt_vals.append(norm_df.loc[curr_seq][tp] - norm_df.loc[mut_seq][tp])
                contexts.append(seq)
                # get the above sequences values and normalize the option containing seq to the mutated, save to opt_vals
            comb_df[''.join(list(opt))] = convert_fc_to_occ(np.array(opt_vals)) # convert opt_vals to occupancy
            # values and save to comb_df

        exp_vals = []
        for i in range(len(comb_df)): # iterate over comb_df row ids
            curr_row = comb_df.iloc[i, :]  # current row
            exp_vals.append(1 - np.prod(1 - np.array(curr_row.values[1:])))  # calculate expected
        comb_df['Expected'] = exp_vals  # add expected values as a  column to the dataframe
        comb_df.rename({''.join(wt_opt): 'Observed', ''.join(comb_opts[1]): '+-', ''.join(comb_opts[2]): '-+'},
                       axis=1, inplace=True) # change wt column to "observed"
        comb_df['Observed-Expected'] = comb_df['Observed'] - comb_df['Expected']  # add a column of the differance
        # between observed and expected
        comb_df['Comb'] = np.repeat(comb_name, len(comb_df))
        comb_df['Context'] = contexts
        all_combs_all_vals = pd.concat([all_combs_all_vals, comb_df])

        mean_obs_comb = comb_df['Observed'].mean() * 100  # mean over Observed
        mean_exp_comb = comb_df['Expected'].mean() * 100  # mean over Expected
        mean_delta_comb = comb_df['Observed-Expected'].mean() * 100  # mean over Observed-Expected
        sem_delta_comb = comb_df['Observed-Expected'].sem() * 100  # calculate sem over Observed-Expected

        # save mean row to mean_comb_df
        curr_comb_row = pd.DataFrame(
            [lib_info['gene'], lib_num, comb_name, min_dist, pos_dist, comb_tfs, mean_obs_comb, mean_exp_comb,
             mean_delta_comb, sem_delta_comb, tf]).T
        curr_comb_row.columns = col_names_mean
        mean_comb_df = pd.concat([mean_comb_df, curr_comb_row], axis=0)

    comb_type = []
    for curr_comb in mean_comb_df['Comb_tfs']:
        split_comb_tfs = curr_comb.split('_')
        tf_comb_pos = [ti for ti, comb_tf in enumerate(split_comb_tfs) if tf[0:3] in comb_tf]
        tf_comb_pos_num = len(tf_comb_pos)
        if tf_comb_pos_num == 2:
            comb_type.append('Self_comb')
        else:
            comb_type.append('non-self_comb')
    mean_comb_df['Comb_type'] = np.reshape(comb_type, (len(mean_comb_df), 1))
    # return mean_comb_df, all_combs_all_vals
    return mean_comb_df


def mot_dist_by_mut_comb(lib_num, curr_comb, pos_by_tf):
    '''This function gets library number, list of 2 str variable position indices (curr_comb)
    and the TFs motifs in each position. The output of this function is the minimal distance
    between the specified tf motifs that contain variable positions in the given ids.'''
    # print(lib_num, curr_comb, pos_by_tf)
    motif_locs = dp.get_motif_info_by_pos(lib_num)  # get df containing tf motif, start and stop position and variable
    # position in motif range.
    tf_mot_pos = []
    for pos_i in range(len(curr_comb)):  # iterate positions indices
        curr_pos = int(curr_comb[pos_i])  # current position
        curr_tf = pos_by_tf[pos_i]  # current tf
        if len(curr_tf) != 3:  # if tf name contains the number
            # curr_tf_split =  re.match(r"([a-z]+)([0-9]+)", curr_tf, re.I).groups()[0].upper()
            # take only the tf letters and upper
            curr_tf_split = curr_tf[0:3].upper()
        else:
            curr_tf_split = curr_tf.upper()  # in case numbers are not included - upper
        pos_subset = motif_locs.query("Pos==@curr_pos")  # subset "curr_tf" by the relevant variable position index
        mot_info = pos_subset.iloc[np.where(pos_subset['TF'].str.contains(curr_tf_split))[0], :]
        # get relevant tf motif info
        # print(tf_mot_pos)
        tf_mot_pos.append([mot_info['Start'].values[0], mot_info['End'].values[0]])
        # store start and stop position in "tf_mot_pos"
    curr_locs = np.reshape(tf_mot_pos, (2, 2))
    # convert the "tf_mot_pos" list to matrix- rows: motifs, columns- start, stop
    first_pos_i = np.argmin(curr_locs[:, 0])  # set the row id of first motif to the motif that starts first
    second_pos_i = np.argmax(curr_locs[:, 0])  # set the row id of second motif to the motif that ends last
    minimal_dist = curr_locs[second_pos_i, 0] - curr_locs[first_pos_i, 1]  # calculate the distance between the
    # two motifs
    if (tf_mot_pos[0] == tf_mot_pos[1]) | (minimal_dist < 0):  # in case that (1) there are 2 positions in same motif,
        # (2) 2 motifs overlap
        return 0  # distance is 0
    else:
        return minimal_dist


def norm_exp_data_tp(tf, lib_num, tp, exp_num, norm_to=None):
    file_name = tf + '_' + str(lib_num) + '_' + str(exp_num) + '.csv'  # find results file ('tf_libnum_exp.csv')
    filt_res = (pd.read_csv(os.path.join(params.RES_PATH, file_name), index_col=0))
    sample_filt, _ = dp.rm_samples(tf, lib_num, exp_num, filt_res)

    if norm_to == 'abs':
        abs_data = dp.get_exp_absolute(tf, lib_num, exp_num)
        abs_filtered, _ = dp.rm_samples(tf, lib_num, exp_num, abs_data)
        log2_norm = sample_filt.div(abs_filtered.values).apply(np.log2)
    else:
        sample_filt_norm = dp.norm_reads(sample_filt)
        log2_norm = dp.res_log2(sample_filt_norm)

    sample_info = dp.get_samp_info(log2_norm)
    norm_tp_0 = dp.norm_to_tp_0(log2_norm, sample_info)

    if norm_to != 'abs':
        norm_tp_0 = norm_non_cut(lib_num, norm_tp_0, norm_to)
    tp_data = norm_tp_0.iloc[:, norm_tp_0.columns.str.contains(tp)]

    return tp_data


def combine_2_exps(tf, lib_num, tp, exps_list, norm_to=None):

    updated_exps_list = []
    for exp_num in exps_list:
        file_name = tf + '_' + str(lib_num) + '_' + str(exp_num) + '.csv'  # find results file ('tf_libnum_exp.csv')
        col_exp = (pd.read_csv(os.path.join(params.RES_PATH, file_name), index_col=0)).columns
        col_tp = [col for col in col_exp if tp in col]
        if len(col_tp) > 0:
            updated_exps_list.append(exp_num)

    if len(updated_exps_list) < 2:
        exp_num = updated_exps_list[0]
        tp_data = norm_exp_data_tp(tf, lib_num, tp, exp_num, norm_to)

    else:
        exps_dfs = []
        for exp_num in updated_exps_list:
            exp_df = norm_exp_data_tp(tf, lib_num, tp, exp_num, norm_to)
            exps_dfs.append(exp_df)
        cols_expA = exps_dfs[0].columns
        rep_i_new = str(int(max([i.split('_')[0] for i in cols_expA]))+1)
        cols_expB = ['_'.join([rep_i_new, col.split('_')[1],  col.split('_')[2]]) for col in exps_dfs[1].columns]
        exps_dfs[1].columns = cols_expB
        comb_data = pd.concat(exps_dfs, axis=1)
        tp_data = comb_data.loc[:, comb_data.columns.str.contains(tp)]
    return tp_data


def get_tf_positions(lib_info, tf):
    '''This function gets the positions of tf motifs and the positions of
    non tf motifs withing a given library.'''
    tf_pos = [i for i, x in enumerate(lib_info['mut_by_tf']) if tf[0:3] in x]
    non_tf_pos = np.setdiff1d(np.array(range(len(lib_info['mut_loc']))), tf_pos)
    return tf_pos, non_tf_pos


def get_min_dist_tf(lib_num, tf, tf1_pos, tf2_pos):
    lib_info = dp.get_lib_info(lib_num)
    comb_tfs = lib_info['mut_by_tf'][tf1_pos] + '_' + lib_info['mut_by_tf'][tf2_pos]  # TFs in combination
    split_comb_tfs = comb_tfs.split('_')

    tf_pos_i = [pos_i for pos_i, tfs in enumerate(comb_tfs.split('_')) if re.findall('\d*\D+', tf)[0] in tfs][0]
    second_tf_pos_i = abs(tf_pos_i - 1)

    comb = [tf1_pos, tf2_pos]
    pos_by_tf = np.array(['xxxx'] * len(comb))
    pos_by_tf[tf_pos_i] = tf
    if '/' in split_comb_tfs[second_tf_pos_i]:
        second_mot_tfs = split_comb_tfs[second_tf_pos_i].split('/')
        min_dist = 1000
        for curr_tf in second_mot_tfs:
            pos_by_tf[second_tf_pos_i] = curr_tf
            curr_dist = mot_dist_by_mut_comb(lib_num, comb, pos_by_tf)
            min_dist = min(min_dist, curr_dist)
    else:
        pos_by_tf[second_tf_pos_i] = split_comb_tfs[second_tf_pos_i]
        min_dist = mot_dist_by_mut_comb(lib_num, comb, pos_by_tf)
    return min_dist


def comb_opts_by_context(tf, lib_num, exp_num, comb, tp):
    lib_info = dp.get_lib_info(lib_num)  # get library information
    wt_var_seq = ''.join(lib_info['wt_at_var_loc'])  # wt seq
    mut_var_seq = get_mut_var(lib_num)  # mutated seq

    norm_df = dp.mean_over_bio_reps(norm_res_data(tf, lib_num, exp_num, norm_to=0))  # get the mean normalized data in
    # fold change values
    seqs = norm_df.index  # library sequences
    all_seqs_split = np.array([list(seq) for seq in seqs])  # library sequences split

    comb_size = len(comb)

    comb_rm_seqs = np.delete(all_seqs_split, np.array(comb), axis=1)  # remove combination from "all_seqs_split"
    flanking_pos_seq_opts = np.unique([''.join(seq) for seq in comb_rm_seqs])  # join all remaining non-comb
    # positions to string sequences
    wt_opt = list(np.array(list(wt_var_seq))[np.array(comb)])  # get the wt option for combination
    mut_opt = list(np.array(list(mut_var_seq))[np.array(comb)])  # get the mutated option for combination
    comb_pos_seq = [[mut_var_seq[pos], wt_var_seq[pos]] for pos in comb]
    # generate a list containing for each position a sublist os the mutated option and the wt option

    #set constant order
    comb_opts = gf.comb_2_arrays(comb_pos_seq[0], comb_pos_seq[1])  # generate all permutations of position option
    # sequence (++, --, +-, -+)
    mut_opt_i = [i for i, opt in enumerate(comb_opts) if list(opt) == mut_opt][0]  # get the id of the mutated
    # option in "comb_opta"
    wt_opt_i = [i for i, opt in enumerate(comb_opts) if list(opt) == wt_opt][0]  # get the id of the wt option
    # in "comb_opts"
    comb_opts_ids = list(range(len(comb_opts)))  # the number of possible combination
    comb_opts_ids.pop(wt_opt_i)  # remove the wt option id from "comb_opts_ids"
    comb_opts_ids.pop(mut_opt_i)  # remove the wt option id from "comb_opts_ids"
    opts_order = [wt_opt_i, 0, 0, mut_opt_i]  # fix position order (++, +-, -+)
    curr_opt = comb_opts[comb_opts_ids[0]]  # get the first combination opts that are not wt or mutated
    if curr_opt[0] == wt_opt[0]:  # if the first position in the option is wt
        opts_order[1] = comb_opts_ids[0]  # this option will be the second item in opts_order
        opts_order[2] = comb_opts_ids[1]  # and the other option will be the third
    else:  # if the first position in the option is mutated
        opts_order[2] = comb_opts_ids[0]  # this option will be the third item in opts_order
        opts_order[1] = comb_opts_ids[1]  # and the other option will be the second
    comb_opts = np.array(comb_opts)[np.array(opts_order)]  # order comb_opts by opts_order


    comb_df = pd.DataFrame()  # generate an empty data frame for the current combination
    for opt in comb_opts:  # iterate over options
        opt_vals = []
        contexts = []
        for seq in flanking_pos_seq_opts:  # iterate over context position sequences
            curr_seq = list(seq)
            for pos in range(comb_size):  # iterate over position ids in combination
                curr_seq.insert(comb[pos], opt[pos])  # insert current option in current context
            curr_seq = ''.join(curr_seq)
            opt_vals.append(norm_df.loc[curr_seq][tp])
            contexts.append(seq)
            # get the above sequences values and normalize the option containing seq to the mutated, save to opt_vals
        comb_df[''.join(list(opt))] = np.array(opt_vals)  # convert opt_vals to occupancy values and save to comb_df
    comb_df.columns = ['++', '+-', '-+', '--']
    comb_df['Context'] = contexts
    return comb_df


if __name__ == '__main__':
    x = 1