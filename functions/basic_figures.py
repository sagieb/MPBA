from functions import params
from functions import data_processing as dp
from functions import lib_table_operations as lto
from functions import general_functions as gf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import OrderedDict
import re


def scatter_mot_mean(res_t, num_full_mot, tf, lib):
    '''This function plots scatter plots of sub-sequences of only variable positions located within tf motifs (and averaged over other positions
    dot color indicates the number of intact tf motifs
    the function inputs include a results table, a vector containing the information of the amount of full tf motif in each sub sequence, the tf library and experiment'''
    curr_lib_info = dp.get_lib_info(lib)  # get a dict of library information
    samp_info = dp.get_samp_info(res_t)  # get a matrix containing for each sample 1:bio_rep 2:time_point 3:tech_Rep
    tps = np.unique(samp_info[:, 1])  # get time points
    num_tps = len(tps)  # get the number of time points
    bio_reps = np.unique(samp_info[:, 0])  # get bio reps
    num_bio_reps = len(bio_reps)  # get the number of bio reps
    num_tech_reps = np.max(samp_info[:, 2])  # get the number of tech reps
    sum_mot_df = pd.DataFrame(num_full_mot)  # turning the vector of amount of full motifs into dataframe
    sum_mot_df.index = res_t.index  # changing the index to go along with the results table
    # we had a few kinds of experiments. in some we had 2/4 tech reps so it easy to compare. in others we had three.
    # in the latter ones we compare each of the three to the two others and hence the if
    if num_tech_reps%2:  # for an odd number of repeats
        if num_tps > 1:  # removing empty rows in figure
            fig, axs = plt.subplots(nrows=num_tps, ncols=num_tech_reps, figsize=(12, 6))
        else:
            fig, axs = plt.subplots(nrows=num_tps, ncols=num_tech_reps, figsize=(12, 3))
        plt.subplots_adjust(hspace=0.3, wspace=1)  # adjusting white space between subplots
        fig.suptitle(tf + ' on ' + curr_lib_info['gene'])  # setting the figure title
        tech_tp_list = list(np.array(range(num_tech_reps)))  # getting a list of technical repeats
        all_pair_combs = [[a, b] for a in tech_tp_list for b in tech_tp_list[a + 1:]] # generating of possible combinations of tech reps
        count = 0  # for getting the wanted ax in axs
        for tp in tps:  # running over time points
            res_tp_filt = res_t.iloc[:, res_t.columns.str.contains(str(tp))]  # subseting the results table taking only columns of the current time point
            for comb in all_pair_combs:  # running over all possible tech rep combinations
                df_for_scat = pd.concat([res_tp_filt.iloc[:, comb[0]], res_tp_filt.iloc[:, comb[1]], sum_mot_df], axis=1)  # building a df for plot
                df_for_scat.columns = [res_tp_filt.columns[comb[0]], res_tp_filt.columns[comb[1]], 'TF mots']  # changing col names for legend
                sns.scatterplot(df_for_scat.iloc[:, 0], df_for_scat.iloc[:, 1], ax=axs.ravel()[count],
                               hue=df_for_scat.iloc[:, 2], s=100, palette='viridis', legend='auto')  # scatter plot!
                sns.move_legend(axs.ravel()[count], "upper left", bbox_to_anchor=(1, 1))  # moving legend to the side
                count += 1
    else:  # for an even number of tech reps
        if num_tps > 1:  # removing empty rows in figure
            fig, axs = plt.subplots(nrows=num_tps, ncols=num_bio_reps, figsize=(9, 12))
        else:
            fig, axs = plt.subplots(nrows=num_tps, ncols=num_tech_reps, figsize=(10, 3))
        plt.subplots_adjust(hspace=0.3, wspace=1)
        fig.suptitle(tf + ' on ' + curr_lib_info['gene'] + ' average over motif combinations')
        count = 0
        for ax in axs.ravel():  #running over axs
            df_for_scat = pd.concat([res_t.iloc[:, count:count+2], sum_mot_df], axis=1)  # taking only current pair tech reps
            df_for_scat.columns = [res_t.columns[count], res_t.columns[count+1], 'TF mots']  # changing col names for legend
            sns.scatterplot(x=df_for_scat.iloc[:, 0], y=df_for_scat.iloc[:, 1],
                            hue=df_for_scat.iloc[:, 2], s=100, palette='viridis', ax=ax, legend='auto')  # scatter plot!
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))  # moving legend to the side
            count += 2


def scatter_tech_reps(res_t, tf, lib, exp):
    '''This function plots scatter plots comparing all technical repeats of tf binding over a given library over all time-points.
    the dot size indicate the number of intact tf motifs in each sequence, the dot color indicates the number of intact wt nuclotides
    in the remaining variable positions (the ones not included in the motifs of the tf).'''
    curr_lib_info = dp.get_lib_info(lib)  # get a dict of library information
    oc_ic = lto.seq_mut_wt_cmp(res_t, curr_lib_info, tf) # for each sequence this function gets the color and size (described above)
    samp_info = dp.get_samp_info(res_t)  # get a matrix containing for each sample 1:bio_rep 2:time_point 3:tech_Rep
    tps = np.unique(samp_info[:, 1])  # get a list of experiment time_points
    num_tps = len(tps)  # get the number of tps
    bio_reps = np.unique(samp_info[:, 0])  # get bio reps
    num_bio_reps = len(bio_reps)  # get the number of bio reps
    num_tech_reps = np.max(samp_info[:, 2])  # get the number of tech reps
    # we had a few kinds of experiments. in some we had 2/4 tech reps so it easy to compare. in others we had three.
    # in the latter ones we compare each of the three to the two others and hence the if
    if num_tech_reps%2:  # for an odd number of tech reps
        if num_tps > 1:  # setting the figure size based on the number of time points
            ncols = num_tech_reps
            fig, axs = plt.subplots(nrows=num_tps, ncols=ncols, figsize=(12, 6))
        else:
            ncols = num_tech_reps
            fig, axs = plt.subplots(nrows=num_tps, ncols=ncols, figsize=(12, 3))
        plt.subplots_adjust(hspace=0.3, wspace=1)  # adjusting white space between sub figures
        fig.suptitle(tf + ' on ' + curr_lib_info['gene'] + '(lib number:' + str(lib) + ' exp:' + str(exp)) # setting the figure title
        tech_tp_list = list(np.array(range(num_tech_reps)))  # getting a list of tech reps
        all_pair_combs = [[a, b] for a in tech_tp_list for b in tech_tp_list[a + 1:]]  # finding all possible combinations of tech reps
        count = 0  # for running over axs
        for tp in tps:  # running over time points
            res_tp_filt = res_t.iloc[:, res_t.columns.str.contains(str(tp))]  # subseting the res table taking only the tech reps that belong to the current time point
            for comb in all_pair_combs:  # plotting all possible pair combinations
                df_for_scat = pd.concat([res_tp_filt.iloc[:, comb[0]], res_tp_filt.iloc[:, comb[1]], oc_ic], axis=1)  # building a dataframe for plotting.
                # the df contains the two comapred tech reps and the data needed for color and size
                df_for_scat.columns = [res_tp_filt.columns[comb[0]], res_tp_filt.columns[comb[1]], 'WT nucs', 'TF mots']  # changing col names for plot
                sns.scatterplot(df_for_scat.iloc[:, 0], df_for_scat.iloc[:, 1], ax=axs.ravel()[count],
                               hue=df_for_scat.iloc[:, 2], size=df_for_scat.iloc[:, 3]/10, palette='viridis', legend='auto', sizes=(20, 200))  # scatter plot!
                # for some reason there's a bug using the size input (it plots it in reverse order) so we divided by 10.
                sns.move_legend(axs.ravel()[count], "upper left", bbox_to_anchor=(1, 1))  # moving the legend to the side
                count += 1
    else:  # for an even number of tech reps
        if num_tps > 1: # setting the figure size based on the number of time points
            ncols = num_bio_reps
            fig, axs = plt.subplots(nrows=num_tps, ncols=ncols, figsize=(9, 12))
        else:
            ncols = num_tech_reps
            fig, axs = plt.subplots(nrows=num_tps, ncols=ncols, figsize=(10, 3))
        plt.subplots_adjust(hspace=0.3, wspace=1)
        fig.suptitle(tf + ' on ' + curr_lib_info['gene'] + '(lib number:' + str(lib) + ' exp:' + str(exp)) # setting the figure title
        count = 0  # for subseting reps
        for ax in axs.ravel():  # running over axs
            df_for_scat = pd.concat([res_t.iloc[:, count:count+2], oc_ic], axis=1) # subseting the compared tech reps from full table and adding color and size data
            df_for_scat.columns = [res_t.columns[count], res_t.columns[count+1], 'WT nucs', 'TF mots'] # changing col names for legened
            sns.scatterplot(x=df_for_scat.iloc[:, 0], y=df_for_scat.iloc[:, 1],
                            hue=df_for_scat.iloc[:, 2], size=df_for_scat.iloc[:, 3]/10, palette='viridis', ax=ax, legend='auto', sizes=(20, 200))  # scatter plot!
            sns.move_legend(ax,"upper left", bbox_to_anchor=(1, 1))  # moving legend to the side
            count+=2
    return fig, num_tps, num_tech_reps

def plot_corr_all_proms(tf):
    '''This function takes a tf and plots the correlation matrix of all samples in all promoters.'''
    dir_files = os.listdir(params.RES_PATH)  #get the list of all csv data files
    tf_files = [file for file in dir_files if tf+'_' in file]  #get only files containing data of the chosen tf. btw  - the "'_'" is the prevent
    # taking Msn2OE (for example) when looking for Msn2
    sqrt_exps = int(np.ceil(np.sqrt(len(tf_files))))
    if sqrt_exps < 2:
        sqrt_exps = 2

    # to chose the number of subplots in the figure we use the square roots of the number of relevant experiments of the tf
    if np.square(sqrt_exps)-sqrt_exps >= len(tf_files):  # because we use ceil in the previous row we're removing empty rows using this if
        fig, axs = plt.subplots(nrows=sqrt_exps-1, ncols=sqrt_exps, figsize=(sqrt_exps * 3, (sqrt_exps - 1) * 4))  #if a row is empty we remove one
    else:
        fig, axs = plt.subplots(nrows=sqrt_exps, ncols=sqrt_exps, figsize=(sqrt_exps*3, sqrt_exps*4))  #if the last row is not empty we keep it
    plt.subplots_adjust(hspace=0.8, wspace=1)
    for i, file in enumerate(tf_files):  #going over all tf files
        _, lib, exp = file.split('.')[0].split('_')  # saving the library number and experiment number
        curr_lib_info = dp.get_lib_info(lib)  # get library information using our get lib info function
        norm_res = dp.norm_reads(dp.filter_results(tf, lib, exp))  # filter right sequences and normalize
        norm_samp_filt, _ = dp.rm_samples(tf, lib, exp, norm_res)
        sns.heatmap(norm_samp_filt.corr(), cmap='viridis', ax=axs.ravel()[i]).set(title=tf+' on ' +
                                                                                        curr_lib_info['gene'] + ' exp:' + exp) # plot correlation heatmap.
    return fig


def scatter_mean(res_t, tf, lib, exp, ax):
    '''This function scatter plots comparing the average of the earliest to the average of the latest time point
    (for a given tf in a given library).
    for cases where there's only one time point (and 3 tech reps) the function compares the first tech rep to an average of 2 & 3.
    the dot size indicate the number of intact tf motifs in each sequence, the dot color indicates the number of intact wt nuclotides
    in the remaining variable positions (the ones not included in the motifs of the tf).'''
    curr_lib_info = dp.get_lib_info(lib) # get a dict of library information
    oc_ic = lto.seq_mut_wt_cmp(res_t,curr_lib_info,tf) # for each sequence this function gets the color and size (described above)
    samp_info = dp.get_samp_info(res_t) # get a matrix containing for each sample 1:bio_rep 2:time_point 3:tech_Rep
    tps = np.unique(samp_info[:,1]) # get a list of experiment time_points
    num_tps = len(tps) # get the number of tps
    bio_reps = np.unique(samp_info[:,0]) # get bio reps
    num_bio_reps = len(bio_reps) # get the number of bio reps
    num_tech_reps = np.max(samp_info[:,2]) # get the number of tech reps
    if num_tps>1: # for most experiments
        early_tp = res_t.iloc[:,res_t.columns.str.contains(str(np.min(tps)))].mean(axis=1) # averaging over the earliest time point
        late_tp = res_t.iloc[:,res_t.columns.str.contains(str(np.max(tps)))].mean(axis=1) # averaging over the latest time poinst
        df_for_scat = pd.concat([early_tp,late_tp,oc_ic],axis=1) # building df adding color and size data
        df_for_scat.columns = [str(np.min(tps)),str(np.max(tps)),'WT nucs','TF mots'] # changing col names for legened
        sns.scatterplot(x=df_for_scat.iloc[:,0],y=df_for_scat.iloc[:,1],
                                hue=df_for_scat.iloc[:,2],size=df_for_scat.iloc[:,3]/10,palette='viridis', legend='auto',ax=ax, sizes=(20,200)) # scatter plot!
        ax.set_title(tf + ' on ' + curr_lib_info['gene'] + ' exp:' + str(exp)) # setting the title
        sns.move_legend(ax,"upper left", bbox_to_anchor=(1, 1)) # moving legend to the side
    elif num_tech_reps*num_bio_reps>3: # for experiments with only 1 tp and 2 bio reps
        df_for_scat = pd.concat([res_t.iloc[:,0:1].mean(axis=1),res_t.iloc[:,-2:].mean(axis=1),oc_ic],axis=1) # averaging over the 1 & 2 time poinst and comparing the result to the 3 & 4
        df_for_scat.columns = [str(np.min(tps)) + '_1&2', str(np.min(tps)) + '_3&4','WT nucs','TF mots'] # changing col names for legened
        sns.scatterplot(x=df_for_scat.iloc[:,0],y=df_for_scat.iloc[:,1],
                                hue=df_for_scat.iloc[:,2],size=df_for_scat.iloc[:,3]/10,palette='viridis', legend='auto',ax=ax, sizes=(20,200)) # scatter plot!
        sns.move_legend(ax,"upper left", bbox_to_anchor=(1, 1)) # moving legend to the side
        ax.set_title(tf + ' on ' + curr_lib_info['gene'] + ' exp:' + str(exp)) # setting the title
    else: # for the experiment where we lost one time point and have only tech reps in one tp
        df_for_scat = pd.concat([res_t.iloc[:,0],res_t.iloc[:,-2:].mean(axis=1),oc_ic],axis=1) # averaging over the 2 & 3 time poinst and comparing the result to the first tp
        df_for_scat.columns = [str(tps[0]) + '_1',str(tps[0]) + '_2&3','WT nucs','TF mots'] # changing col names for legened
        sns.scatterplot(x=df_for_scat.iloc[:,0],y=df_for_scat.iloc[:,1],
                                hue=df_for_scat.iloc[:,2],size=df_for_scat.iloc[:,3]/10,palette='viridis', legend='auto',ax=ax, sizes=(20,200)) # scatter plot!
        sns.move_legend(ax,"upper left", bbox_to_anchor=(1, 1)) # moving legend to the side
        ax.set_title(tf + ' on ' + curr_lib_info['gene'] + ' exp:' + str(exp)) # setting the title


def scatter_mean_mot_mean(res_t,num_mot,tf,lib,exp,ax):
    '''This function generates scatter plots comparing the average of the earliest to the average of the latest time point
    on the sub-sequences composed of positions localized within the motifs of the tf (for a given tf in a given library).
    This function is a bit stupid and the if/elif/else at the end is hard coded and only fits the current structure of out data and experiments'''
    num_mot = pd.Series(num_mot) # turning the number of complete tf motifs in each sub-sequence into a dataframe to be used later in the plot
    num_mot.index = res_t.index # adjusting the index to fit those of the results table
    curr_lib_info = dp.get_lib_info(lib) # get a dict of library information
    samp_info = dp.get_samp_info(res_t) # get a matrix containing for each sample 1:bio_rep 2:time_point 3:tech_Rep
    tps = np.unique(samp_info[:,1]) # get a list of experiment time_points
    num_tps = len(tps) # get the number of tps
    bio_reps = np.unique(samp_info[:,0]) # get bio reps
    num_bio_reps = len(bio_reps) # get the number of bio reps
    num_tech_reps = np.max(samp_info[:,2]) # get the number of tech reps
    if num_tps>1: # for most experiments
        early_tp = res_t.iloc[:,res_t.columns.str.contains(str(np.min(tps)))].mean(axis=1) # averaging over the earliest time point
        late_tp = res_t.iloc[:,res_t.columns.str.contains(str(np.max(tps)))].mean(axis=1) # averaging over the latest time poinst
        df_for_scat = pd.concat([early_tp,late_tp,num_mot],axis=1) # building df adding color and size data
        df_for_scat.columns = [str(np.min(tps)),str(np.max(tps)),'TF mots'] # changing col names for legened
        sns.scatterplot(x=df_for_scat.iloc[:,0],y=df_for_scat.iloc[:,1],
                                hue=num_mot,palette='viridis', legend='auto',ax=ax, sizes=(20,200)) # scatter plot!
        ax.set_title(tf + ' on ' + curr_lib_info['gene'] + ' exp:' + str(exp)) # setting the title
        sns.move_legend(ax,"upper left", bbox_to_anchor=(1, 1)) # moving legend to the side
    elif num_tech_reps*num_bio_reps>3: # for experiments with only 1 tp and 2 bio reps
        df_for_scat = pd.concat([res_t.iloc[:,0:1].mean(axis=1),res_t.iloc[:,-2:].mean(axis=1),num_mot],axis=1) #  averaging over the 1 & 2 time poinst and comparing the result to the 3 & 4
        df_for_scat.columns = [str(tps[0]) + '_1&2', str(tps[0]) + '_3&4','TF mots'] # changing col names for legened
        sns.scatterplot(x=df_for_scat.iloc[:,0],y=df_for_scat.iloc[:,1],
                                hue=num_mot,palette='viridis', legend='auto',ax=ax, sizes=(20,200)) # scatter plot!
        sns.move_legend(ax,"upper left", bbox_to_anchor=(1, 1)) # moving legend to the side
        ax.set_title(tf + ' on ' + curr_lib_info['gene'] + ' exp:' + str(exp)) # setting the title
    else: # for the experiment where we lost one time point and have only tech reps in one tp
        df_for_scat = pd.concat([res_t.iloc[:,0],res_t.iloc[:,-2:].mean(axis=1),num_mot],axis=1) # averaging over the 2 & 3 time poinst and comparing the result to the first tp
        df_for_scat.columns = [str(tps[0]) + '_1',str(tps[0]) + '_2&3','TF mots'] # changing col names for legened
        sns.scatterplot(x=df_for_scat.iloc[:,0],y=df_for_scat.iloc[:,1],
                                hue=num_mot,palette='viridis', legend='auto',ax=ax, sizes=(20,200)) # scatter plot!
        sns.move_legend(ax,"upper left", bbox_to_anchor=(1, 1)) # moving legend to the side
        ax.set_title(tf + ' on ' + curr_lib_info['gene'] + ' exp:' + str(exp)) # setting the title


def scatter_mean_tf_proms(tf, mean_over_mot=None):
    '''This function takes a tf and scatterplots the change of each sequence - if 'mean_over_mot' is empty then the function will use scatter_mean
    to plot the mean of the earliest time point on one axes and the latest on the other axis for each promoter. if mean_over_mot=True then the function will use 'scatter_mean_mot_mean'
    and plot the same only for averaged sub-sequences of variable positions localized within tf motifs'''
    # for all other functions in the script
    dir_files = os.listdir(params.RES_PATH) #get the list of all csv data files
    tf_files = [file for file in dir_files if tf + '_' in file] #get only files containing data of the chosen tf. btw  - the "'_'" is the prevent
    # taking Msn2OE (for example) when looking for Msn2
    sqrt_exps = int(np.ceil(np.sqrt(len(tf_files)))) # to chose the number of subplots in the figure we use the square roots of the number of relevant experiments of the tf
    if sqrt_exps < 2:
        sqrt_exps =2
    if np.square(sqrt_exps)-sqrt_exps>=len(tf_files): # because we use ceil in the previous row we're removing empty rows using this if
        fig, axs = plt.subplots(nrows=sqrt_exps-1, ncols=sqrt_exps, figsize=((sqrt_exps)*4,(sqrt_exps-1)*4)) #if a row is empty we remove one
    else:
        fig, axs = plt.subplots(nrows=sqrt_exps, ncols=sqrt_exps, figsize=(sqrt_exps*4,sqrt_exps*4)) #if the last row is not empty we keep it
    plt.subplots_adjust(hspace=0.8, wspace=1.8)
    for i, file in enumerate(tf_files): #going over all tf files
        _, lib, exp = file.split('.')[0].split('_') # saving the library number and experiment number
        norm_res = dp.norm_reads(dp.filter_results(tf, lib, exp)) # filter right sequences and normalize
        _, sample_filt_norm = dp.rm_samples(tf, lib, exp, norm_res)
        log2_norm = dp.res_log2(sample_filt_norm)
        sample_info_with_tp0 = dp.get_samp_info(log2_norm)  # getting a matrix with experiment information
        res_t = dp.norm_to_tp_0(log2_norm, sample_info_with_tp0)  # generating a tp 0 normalized table
        norm_uncut = lto.norm_non_cut(res_t, 3)
        if mean_over_mot:
            norm_uncut,num_mot = lto.get_mot_mean_table(norm_uncut, lib, tf) #get averaged sub-sequence table and the number of full motif in each
            scatter_mean_mot_mean(norm_uncut, num_mot, tf, lib, exp, axs.ravel()[i])  # scattering using the 'scatter_mean_mot_mean' giving it the right ax to plot on.
        else:
            scatter_mean(norm_uncut, tf, lib, exp, axs.ravel()[i]) # scattering using the scatter mean function giving it the right ax to plot on.
    return fig

def plot_abs_on_prom(lib_num):
    '''This function plots the time point averaged log2 fold change of all tfs compared to the absolute measurement on a chosen library.'''
    # abs_exps = [12,13,14,16] #experiments with absolute measurements
    abs_exps = [13, 14, 16] # removed exp 12 for now because its shitty
    lib_files = dp.find_lib_files(lib_num)  # get all files of the chosen library
    lib_info = dp.get_lib_info(lib_num)  # get library info
    fin_df = pd.DataFrame()  # opening an empty dataFrame
    for file_name in lib_files:  # running over library files
        tf,_,exp = file_name.split('.')[0].split('_')  # getting the current tf and experiment
        if int(exp) in abs_exps:  # not all experiments have absolute measurements. go in only if the exp was chosen and listed in 'abs_exps'
            res_table = dp.filter_results(tf, lib_num, exp)  # filter relevant sequences without normalization!
            abs_tf_lib_exp = dp.get_exp_absolute(tf, lib_num, exp)  # get absolute reads in each sample
            res_norm_abs = res_table.sum()/abs_tf_lib_exp  # normalize the number of reads in each sample to the number of 'absolute' reads
            sample_info_with0 = dp.get_samp_info(res_norm_abs)  #get sample information (bio_rep,tp,tech_rep)
            abs_ratio_table = dp.norm_to_tp_0(res_norm_abs, sample_info_with0, space='linear') # normalize the absolute normalized res table to time zero
            sample_info_no0 = dp.get_samp_info(abs_ratio_table)  #get sample information (bio_rep,tp,tech_rep) without time zero
            tps = np.unique(sample_info_no0[:, 1])  # get activated time points
            # We now want to average all repeats of a given time point:
            tf_list = [tf]*len(tps)  # generate a list with repeating tf name corresponding to length of 'tps'
            exp_list = [exp]*len(tps)  # same for experiment
            col_names = [m+'_'+str(n)+'_'+str(l) for m, n, l in zip(tf_list, tps, exp_list)] # generating column names for time point averaged table
            for i, tp in enumerate(tps):  # running over activated time points
                curr_df = pd.DataFrame(abs_ratio_table.iloc[:, abs_ratio_table.columns.str.contains(str(tp))].apply(np.log2).mean(axis=1), columns=[col_names[i]])
                # a = np.mean(np.log2(abs_ratio_table.iloc[:, abs_ratio_table.columns.str.contains(str(tp))]).values)
                # print(a)

                # curr_df = pd.DataFrame(a,columns=[col_names[i]])
                # getting the averaged values over all samples from the same time point ^. we calculate 1-values to get the fraction cut!
                # The steps for getting here were: 1. Divide the sum of recognized promoter sequence reads recieved for a given sample by the number
                # of absolute reads. 2. Normalizing this results by time zero. 3. log2 transformation 4. Averaging over all samples of a time point and calcultaing
                # 1 - answer to get the fraction cut.
                fin_df = pd.concat([fin_df, curr_df], axis=1)  # concatanating tf/time point into dataframe
    fig = plt.figure(figsize=(10,8)) # opening new figure
    plt.rcParams.update({'font.size': 22})  # changing font size for the figure
    (fin_df).iloc[0,:].plot.bar(color='#9776E3')  # plotting a bar plot to show the % cut for each tf in each time_point. multiplying by 100 to get to %.
    plt.title('Absolute binding of ' + lib_info['gene'] +' (lib - ' + str(lib_num) + ')')  # setting the figure title
    plt.ylabel('log2(fold change)') # setting ylabel
    return fig


def plot_abs_for_tf(tf):
    '''This function plots the time point averaged log2 fold change of all tf promoters compared to the absolute measurement.'''
    abs_exps = [12,13,14,16] #experiments with absolute measurements
    # abs_exps = [13,14,16] # removed exp 12 for now because its shitty
    tf_files = dp.find_tf_files(tf) # get all files of the chosen tf
    # lib_info = dp.get_lib_info(lib) # get library info
    fin_df = pd.DataFrame() # opening an empty dataFrame
    for file_name in tf_files: # running over tf files
        _,lib_num,exp = file_name.split('.')[0].split('_') # getting the current library and experiment
        if int(exp) in abs_exps: # not all experiments have absolute measurements. go in only if the exp was chosen and listed in 'abs_exps'
            lib_info = dp.get_lib_info(lib_num) # get library info
            res_table = dp.filter_results(tf,lib_num,exp) # filter relevant sequences without normalization!
            abs_tf_lib_exp = dp.get_exp_absolute(tf,lib_num,exp) # get absolute reads in each sample
            res_norm_abs = res_table.sum()/abs_tf_lib_exp # normalize the number of reads in each sample to the number of 'absolute' reads
            sample_info_with0 = dp.get_samp_info(res_norm_abs) #get sample information (bio_rep,tp,tech_rep)
            abs_ratio_table = dp.norm_to_tp_0(res_norm_abs,sample_info_with0,space='linear') # normalize the absolute normalized res table to time zero
            sample_info_no0 = dp.get_samp_info(abs_ratio_table) #get sample information (bio_rep,tp,tech_rep) without time zero
            tps = np.unique(sample_info_no0[:,1]) # get activated time points
            # We now want to average all repeats of a given time point:
            tf_list = [tf]*len(tps) # generate a list with repeating tf name corresponding to length of 'tps'
            exp_list = [exp]*len(tps) # same for experiment
            col_names = [m+'_'+str(n)+'_'+str(l)+'_'+lib_info['gene'] +' (lib - ' + str(lib_num) + ')' for m,n,l in zip(tf_list,tps,exp_list)] # generating column names for
            # time point averaged table with promoter name
            for i,tp in enumerate(tps): # running over activated time points
                # curr_df = pd.DataFrame(np.log2(abs_ratio_table.iloc[:,abs_ratio_table.columns.str.contains(str(tp))].mean(axis=1)),columns=[col_names[i]])
                curr_df = pd.DataFrame(abs_ratio_table.iloc[:, abs_ratio_table.columns.str.contains(str(tp))].apply(np.log2).mean(axis=1), columns=[col_names[i]])
                # building curr_df with log2 on the averaged time point.
                # getting the averaged values over all samples from the same time point ^. we calculate 1-values to get the fraction cut!
                # The steps for getting here were: 1. Divide the sum of recognized promoter sequence reads recieved for a given sample by the number
                # of absolute reads. 2. Normalizing this results by time zero. 3. log2 transformation 4. Averaging over all samples of a time point and calcultaing
                # 1 - answer to get the fraction cut.
                fin_df = pd.concat([fin_df,curr_df],axis=1) # concatanating tf/time point into dataframe
    fig = plt.figure(figsize=(10,8)) # opening new figure
    plt.rcParams.update({'font.size': 12}) # changing font size for the figure
    (fin_df.iloc[0,:]).plot.bar(color='#9776E3') # plotting a bar plot to show the % cut for each tf in each time_point. multiplying by 100 to get to %.
    plt.title('Absolute binding of ' + tf) # setting the figure title
    plt.ylabel('log2(fold change)') # setting ylabel
    plt.tight_layout()
    return fig


def scatter_mean_tfs_on_lib(lib_num, mean_over_mot=None):
    '''This function takes a library and scatterplots the change of each sequence - averaging the earliest time point on one axes and the latest on the other axis for each tf.'''
    lib_files = dp.find_lib_files(lib_num) # get all result files from the chosen lib
    sqrt_exps = int(np.ceil(np.sqrt(len(lib_files))))  #get sqrt exps for number of sub figures
    if sqrt_exps < 2:
        sqrt_exps = 2
    plt.rcParams.update({'font.size': 10})
    if np.square(sqrt_exps) - sqrt_exps >= len(lib_files):  # because we use ceil in the previous row we're removing empty rows using this if
        fig, axs = plt.subplots(nrows=sqrt_exps-1, ncols=sqrt_exps, figsize=((sqrt_exps)*4,(sqrt_exps-1)*4)) #if a row is empty we remove one
    else:
        fig, axs = plt.subplots(nrows=sqrt_exps, ncols=sqrt_exps, figsize=(sqrt_exps*4, sqrt_exps*4)) #if the last row is not empty we keep it
    plt.subplots_adjust(hspace=0.8, wspace=1.8)
    for i,file in enumerate(lib_files): #going over all library files
        tf, _, exp = file.split('.')[0].split('_') # saving the tf and experiment number
        norm_res = dp.norm_reads(dp.filter_results(tf, lib_num, exp)) # filter right sequences and normalize
        _, sample_filt_norm = dp.rm_samples(tf, lib_num, exp, norm_res)
        log2_norm = dp.res_log2(sample_filt_norm)
        sample_info_with_tp0 = dp.get_samp_info(log2_norm)  # getting a matrix with experiment information
        res_t = dp.norm_to_tp_0(log2_norm, sample_info_with_tp0)  # generating a tp 0 normalized table
        norm_uncut = lto.norm_non_cut(res_t, 3)
        if mean_over_mot:
            norm_uncut,num_mot = lto.get_mot_mean_table(norm_uncut, lib_num, tf) #get averaged sub-sequence table and the number of full motif in each
            scatter_mean_mot_mean(norm_uncut, num_mot, tf, lib_num, exp, axs.ravel()[i])  # scattering using the 'scatter_mean_mot_mean' giving it the right ax to plot on.
        else:
            scatter_mean(norm_uncut, tf, lib_num, exp, axs.ravel()[i]) # scattering using the scatter mean function giving it the right ax to plot on.
    return fig

def mean_pos_change(tf, lib_num, exp_num, ax):
    norm_res = dp.norm_reads(dp.filter_results(tf, lib_num, exp_num))  # filter right sequences and normalize
    sample_filt_norm,_ = dp.rm_samples(tf, lib_num, exp_num, norm_res)
    sample_info = dp.get_samp_info(sample_filt_norm)
    mean_tech = dp.mean_bio_rep_tps(sample_filt_norm , sample_info)
    log2_norm = dp.res_log2(mean_tech)
    # mean_repeats = dp.mean_bio_rep_tps(log2_norm, sample_info)
    norm_tp_0 = dp.norm_to_tp_0(log2_norm, sample_info)
    # norm_tp_0 = lto.norm_non_cut(norm_tp_0,3,space=None)
    colors_dict = {'45': 'orange', '60': 'peru', '90': 'sienna', '180': 'maroon'}
    edges = 5
    lib_info_dict = dp.get_lib_info(lib_num)
    mut_loc = lib_info_dict['mut_loc']
    relevant_seqs = list(norm_tp_0.index)
    wt_var_seq = ''.join(lib_info_dict['wt_at_var_loc'])
    col_tps = [x.split('_')[1] for x in norm_tp_0.columns.values]
    tps = list(OrderedDict.fromkeys(col_tps))
    colors = [colors_dict[tp] for tp in tps]
    fig = plt.figure()
    fig.set_size_inches(25, 15)

    for i, curr_opts in enumerate(wt_var_seq):
        relevant_seqs_ids = [j for j, seq in enumerate(relevant_seqs) if seq[i] == curr_opts]
        pos_change = lto.mean_vals_by_id(norm_tp_0, relevant_seqs_ids)
        ax.scatter(np.repeat(mut_loc[i], len(pos_change)), np.array(pos_change), c=colors, edgecolors='k', s=60,
                   label='_nolegend_')
        ax.scatter(mut_loc[i], np.mean(pos_change), color='k', marker='_', s=200, label='_nolegend_')
        only_pos_id = lto.get_only_pos_intact([i], wt_var_seq, norm_tp_0)
        only_pos_change = lto.mean_vals_by_id(norm_tp_0, only_pos_id)
        ax.scatter(np.repeat(mut_loc[i], len(pos_change)), np.array(only_pos_change), c=colors, marker='x', s=20,
                   label='_nolegend_')

    ax.set_xlim(0, params.LIB_SIZE)
    ax.set_ylabel('Log2 FC')
    ax.set_xlabel('Position')
    ax.axhline(y=0, color='k', ls='--', label='_nolegend_')
    ax.set_xlim(min(mut_loc) - edges, max(mut_loc) + edges)

    pair_combs = gf.get_pair_combs(list(range(len(wt_var_seq))), 2)
    neighbor_pairs = [pair for pair in pair_combs if pair[1] - pair[0] == 1]
    for pair in neighbor_pairs:
        pos1_mean, pos2_mean, both_mean = lto.mean_pair_opts(pair, wt_var_seq, norm_tp_0)
        pos_1_loc = mut_loc[pair[0]]
        pos_2_loc = mut_loc[pair[1]]
        pos_dist = np.abs(pos_2_loc - pos_1_loc)
        ax.scatter(np.repeat(pos_1_loc + pos_dist / 3, len(pos1_mean)), np.array(pos1_mean), c=colors, marker='^',
                   edgecolors='k', s=20, label='_nolegend_')
        ax.scatter(np.repeat(pos_2_loc - pos_dist / 3, len(pos2_mean)), np.array(pos2_mean), c=colors, marker='^',
                   edgecolors='k', s=20, label='_nolegend_')
        ax.scatter(np.repeat(pos_1_loc + pos_dist / 2, len(both_mean)), np.array(both_mean), c=colors, marker='^',
                   edgecolors='k', s=20, label='_nolegend_')

    best_pair = ''
    best_pos1_mean = []
    best_pos2_mean = []
    best_both_mean = [np.inf]
    for pair in pair_combs:
        pos1_mean, pos2_mean, both_mean = lto.mean_pair_opts(pair, wt_var_seq, norm_tp_0)
        if min(both_mean) < min(best_both_mean):
            best_pair = pair
            best_pos1_mean = pos1_mean
            best_pos2_mean = pos2_mean
            best_both_mean = both_mean
    pos_1_loc = mut_loc[best_pair[0]]
    pos_2_loc = mut_loc[best_pair[1]]
    pos_dist = np.abs(pos_2_loc - pos_1_loc)
    x_vals = [pos_1_loc, (pos_1_loc + pos_dist / 2), pos_2_loc]
    for ti in range(len(tps)):
        ax.plot(x_vals, [best_pos1_mean[ti], best_both_mean[ti], best_pos2_mean[ti]], c=colors[ti], label=tps[ti])

    ax.set_ylim(-np.max(np.abs(ax.get_ylim())), np.max(np.abs(ax.get_ylim())))
    legender = []
    for i, curr_mut in enumerate(mut_loc):
        if re.split(r'(\d+)', tf)[0] in lib_info_dict['mut_by_tf'][i]:
            ax.scatter(curr_mut, np.max(ax.get_ylim()) * 0.9, edgecolors='k', s=80, linewidth=2)
        else:
            ax.scatter(curr_mut, np.max(ax.get_ylim()) * 0.9, s=80)

        legender.append(lib_info_dict['mut_by_tf'][i])
    ax.legend(tps + legender, loc='upper left', bbox_to_anchor=(1, 1.1), prop={'size': 7.3})
    ax.set_title('TF: ' + tf + ', Promoter: ' + lib_info_dict['gene'] + ', Experiment: ' + str(exp_num))
    plt.tight_layout()
    return fig


def mean_pos_change_iterate_files(files):
    total_subplot_in_slide = 12
    col_num = 2
    row_num = int(total_subplot_in_slide / 2)
    # fig_num = np.ceil(len(files) / total_subplot_in_slide)

    files_groups = []
    for i in np.arange(0, len(files), total_subplot_in_slide):
        files_groups.append(files[i:i + total_subplot_in_slide])

    figs = []
    for group in files_groups:
        fig, axes = plt.subplots(nrows=row_num, ncols=col_num, figsize=(10 * col_num, 2 * row_num))
        plt.subplots_adjust(hspace=1, wspace=0.5)
        for i, ax in enumerate(axes.ravel()):
            if i < len(group):
                file_name = group[i].split('.csv')[0]
                split_file_name = file_name.split('_')
                tf = split_file_name[0]
                lib_num = split_file_name[1]
                # curr_prom = dp.get_lib_info(lib_num)['gene']
                exp_num = split_file_name[-1]
                mean_pos_change(tf, lib_num, exp_num, ax)
        plt.close()
        figs.append(fig)
    return figs


def corr_by_lib_plot(lib_num):
    lib_files = dp.find_lib_files(lib_num)
    lib_info_dict = dp.get_lib_info(lib_num)
    lib_tfs_signal_norm = pd.DataFrame()
    lib_tfs_signal = pd.DataFrame()
    for file_name in lib_files:
        file_name = file_name.split('.csv')[0]
        split_file_name = file_name.split('_')
        tf = split_file_name[0]
        exp_num = split_file_name[-1]
        norm_res = dp.norm_reads(dp.filter_results(tf, lib_num, exp_num))  # filter right sequences and normalize
        sample_filt_norm, _ = dp.rm_samples(tf, lib_num, exp_num, norm_res)
        sample_info = dp.get_samp_info(sample_filt_norm)
        mean_repeats = dp.mean_bio_rep_tps(sample_filt_norm, sample_info)
        log2_norm = dp.res_log2(mean_repeats)
        norm_tp_0 = dp.norm_to_tp_0(log2_norm, sample_info).sort_index()
        new_col_names = [tf + '_' + col for col in norm_tp_0.columns]
        norm_tp_0.columns = new_col_names
        lib_tfs_signal_norm = pd.concat([lib_tfs_signal_norm, norm_tp_0], axis=1)
        lib_tfs_signal
        new_col_names = [tf + '_' + col for col in log2_norm.columns]
        log2_norm.columns = new_col_names
        lib_tfs_signal = pd.concat([lib_tfs_signal, log2_norm], axis=1)

    fig, axes = plt.subplots(figsize=(12, 6), nrows=1, ncols=2)
    sns.heatmap(lib_tfs_signal_norm.corr(), cmap='viridis', ax=axes[0])
    axes[0].set_title('Time point 0 normalized correlation')
    sns.heatmap(lib_tfs_signal.corr(), cmap='viridis', ax=axes[1])
    axes[1].set_title('Log2 signal correlation')
    plt.suptitle(lib_info_dict['gene'])
    return fig

def plot_heatmap(data, title, col_map, ax=None):
    '''This function plot a heatmap based on the following parameters it gets:
    1. data (matrix/data frame)  2.title and color map  3.and optionally an ax.'''
    sns.heatmap(data, cmap=col_map, annot=True, fmt=".2f",  annot_kws={"size": 18.5 / (np.sqrt(len(data))*2)}, ax=ax)
    ax.set_title(title)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    return ax


def plot_pos_change_mat(lib_num, comb_df, tp, lower_mat_col, upper_mat_col, title, ax):
    '''This function get the comb_df generate by the function lto.generate_pos_change_df, and by specifing the relevant columns
    of the data frame it plots an asimetric data frame showing the values pairwize (posXpos) by the input time point'''
    curr_tp_df = comb_df.loc[comb_df['Time_point']==tp] # narrow the comb_df to the relevant time point data
    lower_mat = np.transpose(pd.pivot_table(curr_tp_df, values=lower_mat_col, index=['Pos1'], columns=['Pos2'])) # generate the lower matrix part
    upper_mat = pd.pivot_table(curr_tp_df, values=upper_mat_col, index=['Pos1'], columns=['Pos2']) # generate the upper matrix part
    mat_to_plot = np.where(np.isnan(lower_mat) & np.isnan(upper_mat),np.nan,np.nansum(np.stack((lower_mat,upper_mat)), axis=0)) # combine matrix parts
    lib_info_dict = dp.get_lib_info(lib_num) # get lib info
    df_to_show = pd.DataFrame(mat_to_plot, index=lib_info_dict['mut_by_tf'], columns=lib_info_dict['mut_by_tf'])
    # generate a data frame based on the matrix, with the relevant TF names as index and columns
    col_map = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True) # generate colormap
    col_map.colors = np.flipud(col_map.colors)
    ax = plot_heatmap(df_to_show, title, col_map, ax)  # plot as heat map
    return ax


def single_pos_change_heatmaps(tf, lib_num, exp_num):
    col_map = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True)
    col_map.colors = np.flipud(col_map.colors)

    mean_change_by_pos = lto.mean_pos_change_calc(tf, lib_num, exp_num)
    only_intact_pos_change = lto.only_pos_intact_change_calc(tf, lib_num, exp_num)
    only_mut_pos_change = lto.only_pos_mut_change_calc(tf, lib_num, exp_num)
    one_pos_mut_rest_mean = lto.one_pos_mut_rest_mean_clac(tf, lib_num, exp_num)

    fig, axes = plt.subplots(nrows=4, ncols=1)
    fig.set_size_inches(4, 12)
    plot_heatmap(mean_change_by_pos, 'Mean change by pos', col_map, axes[0])
    plot_heatmap(only_intact_pos_change, 'Only pos wt change', col_map, axes[1])
    plot_heatmap(only_mut_pos_change, 'Only pos mutated change', col_map, axes[2])
    plot_heatmap(one_pos_mut_rest_mean, 'One pos mutated mean over rest', col_map, axes[3])
    fig.tight_layout(pad=0.3)
    return fig


def pair_pos_change_heatmaps(tf, lib_num, exp_num):
    col_map = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True)
    col_map.colors = np.flipud(col_map.colors)
    lib_info_dict = dp.get_lib_info(lib_num)
    comb_df = lto.generate_pos_change_df(tf, lib_num, exp_num)
    tps = comb_df['Time_point'].unique()

    fig, axes = plt.subplots(nrows=3, ncols=len(tps))
    fig.set_size_inches(12, 14)
    i = 0
    j = 0
    for tp in tps:
        if len(tps) == 1:
            axes_locs = [[i], [i+1], [i+2]]
            title_pos_only = 'Only one pos intact\n' + tf + ' - ' + lib_info_dict['gene'] + ', time point: ' + tp
            plot_pos_change_mat(lib_num, comb_df, tp, 'Pos1_only_mean', 'Pos2_only_mean', title_pos_only, axes[i])
            title_both_mean = 'Lower: both wt mean, Upper: both mutated mean \n' + tf + ' - ' + lib_info_dict[
                'gene'] + ', time point: ' + tp
            plot_pos_change_mat(lib_num, comb_df, tp, 'Both_wt_mean', 'Both_mut_mean', title_both_mean, axes[i +1])
            title_both_seq = 'Lower: both wt seq, Upper: both mutated seq\n' + tf + ' - ' + lib_info_dict[
                'gene'] + ', time point: ' + tp
            plot_pos_change_mat(lib_num, comb_df, tp, 'Both_wt', 'Both_mut', title_both_seq, axes[i + 2])
        else:
            title_pos_only = 'Only one pos intact\n' + tf + ' - ' + lib_info_dict['gene'] + ', time point: ' + tp
            plot_pos_change_mat(lib_num, comb_df, tp, 'Pos1_only_mean', 'Pos2_only_mean', title_pos_only, axes[i, j])
            title_both_mean = 'Lower: both wt mean, Upper: both mutated mean \n' + tf + ' - ' + lib_info_dict[
                'gene'] + ', time point: ' + tp
            plot_pos_change_mat(lib_num, comb_df, tp, 'Both_wt_mean', 'Both_mut_mean', title_both_mean, axes[i + 1, j])
            title_both_seq = 'Lower: both wt seq, Upper: both mutated seq\n' + tf + ' - ' + lib_info_dict[
                'gene'] + ', time point: ' + tp
            plot_pos_change_mat(lib_num, comb_df, tp, 'Both_wt', 'Both_mut', title_both_seq, axes[i + 2, j])
            j += 1

    fig.tight_layout(pad=0.3)
    return fig

def plot_abs_mostleast_bytf(tf, num_seqs):
    '''This function plots the mean change of the top `num_seqs` sequence variants over all promoters in all timepoints for a chosen tf'''
    files = list(np.sort(dp.find_tf_files(tf))) # get tf files
    wanted_cols = ['Most','Least','promoter','timepoint'] # df column names
    fin_df = pd.DataFrame(columns=wanted_cols) # empty final df
    for file in list(np.sort(files)): # run over all tf files
        _,lib_num,exp = file.split('.')[0].split('_') # get curr file library number
        if int(exp) not in [12,18,19]: # skip experiments with no absolute measurements
            lib_info = dp.get_lib_info(lib_num) # get curr library info
            curr_norm_abs = dp.norm_var_to_abs(tf,lib_num,exp) # get absolte normalized dataframe
            most_cut,least_cut = lto.get_extreme_cut(tf,lib_num,exp,num_seqs) # get best and worst cut sequences index
            curr_norm_abs.loc[most_cut].mean() # get mean change of most cut sequence (change compared to absolute)
            curr_norm_abs.loc[least_cut].mean() # get mean change of least cut sequence (change compared to absolute)
            # generating 'promoter' column conataining the promoter name and library number of each sample for each:
            curr_prom_col = pd.Series([(lib_info['gene']+':'+str(lib_num))]*len(curr_norm_abs.loc[most_cut].mean()),index=curr_norm_abs.loc[most_cut].mean().index)
            # generating 'timepoints' column conataining the tp of each sample for each:
            curr_timepoints = pd.Series([int(i.split('_')[1]) for i in curr_norm_abs.loc[most_cut].mean().index.values],index=curr_norm_abs.loc[most_cut].mean().index)
            # generating curr experiment dataframe and calculating the absolute change in % of most cut and least cut
            curr_df = pd.concat([(2**curr_norm_abs.loc[most_cut].mean()-1)*100
                                 ,(2**curr_norm_abs.loc[least_cut].mean()-1)*100,pd.Series(curr_prom_col),curr_timepoints],axis=1)
            curr_df.columns = wanted_cols # setting column names to fit the combined final data frame
            fin_df = pd.concat([fin_df,curr_df]) # concatanating curr experiment df and final df
            if 'Msn2' in tf: # this if removes data of hap4 promoter from Msn2 strains. done because it is cut way less than the absolute
                rm_id = [i for i,x in enumerate(fin_df['promoter']) if 'Hap4' not in x]
                fin_df = fin_df.iloc[rm_id,:]
    if len(fin_df) != 0:
        fin_df['Most-Least'] = fin_df['Most']-fin_df['Least'] #new most minus least column
        sns.set(font_scale=1)
        fig,axs = plt.subplots(1,2,figsize=(16, 8)) #open figure
        plt.suptitle(tf)
        # left plot
        curr_ax = 0
        # scattering results
        sns.scatterplot(data=fin_df,x='Most',y='Least',hue='promoter',size='timepoint',sizes=(50,200),palette='Paired',ax=axs.ravel()[curr_ax])
        # adding x=0 y=0 and x=y lines
        axs.ravel()[curr_ax].axhline(0,c='k')
        axs.ravel()[curr_ax].axvline(0,c='k')
        axs.ravel()[curr_ax].axline((0, 0), slope=1,c='k')
        # setting labels
        axs.ravel()[curr_ax].set_xlabel('Most cut seqs change %', fontsize=18)
        axs.ravel()[curr_ax].set_ylabel('Least cut seqs change %', fontsize=18)
        # removing legend from left plot
        axs.ravel()[curr_ax].legend([],[], frameon=False)
        # right plot (y axis showing most-least)
        curr_ax = 1
        sns.scatterplot(data=fin_df,x='Most',y='Most-Least',hue='promoter',size='timepoint',sizes=(50,200),palette='Paired',ax=axs.ravel()[curr_ax])
        axs.ravel()[curr_ax].axhline(0,c='k')
        axs.ravel()[curr_ax].axvline(0,c='k')
        axs.ravel()[curr_ax].axline((0, 0), slope=1,c='k')
        axs.ravel()[curr_ax].set_xlabel('Most cut seqs change %', fontsize=18)
        axs.ravel()[curr_ax].set_ylabel('Most-Least %', fontsize=18)
        sns.move_legend(axs.ravel()[curr_ax], "upper left", bbox_to_anchor=(1, 1))
        return fig, fin_df


def plot_abs_mostleast_bylib(lib_num,num_seqs):
    '''This function plots the mean change of the top `num_seqs` sequence variants over all tfs in all timepoints for a chosen library'''
    files = list(np.sort(dp.find_tf_files(lib_num))) # getting all files of the chosen library
    wanted_cols = ['Most','Least','tf','timepoint'] # df column names
    fin_df = pd.DataFrame(columns=wanted_cols)  # empty final df
    for file in list(np.sort(files)): # run over all library files
        tf,lib_num,exp = file.split('.')[0].split('_') # get curr file exp info
        lib_info = dp.get_lib_info(lib_num)  # get curr library info
        if (int(exp) not in [12,18,19]) & ('Context' not in tf): #skipping matan's context experiments and experiments with no absolute measurement
            curr_norm_abs = dp.norm_var_to_abs(tf,lib_num,exp) # get absolte normalized dataframe
            most_cut, least_cut = lto.get_extreme_cut(tf,lib_num,exp,num_seqs) # get best and worst cut sequences index
            curr_norm_abs.loc[most_cut].mean() # get mean change of most cut sequence (change compared to absolute)
            curr_norm_abs.loc[least_cut].mean() # get mean change of least cut sequence (change compared to absolute)
            # generating tf column conataining the tf name of each sample for each:
            curr_tf_col = pd.Series([(tf)]*len(curr_norm_abs.loc[most_cut].mean()),index=curr_norm_abs.loc[most_cut].mean().index)
            # generating 'timepoints' column conataining the tp of each sample for each:
            curr_timepoints = pd.Series([int(i.split('_')[1]) for i in curr_norm_abs.loc[most_cut].mean().index.values],index=curr_norm_abs.loc[most_cut].mean().index)
            # generating curr experiment dataframe and calculating the absolute change in % of most cut and least cut
            curr_df = pd.concat([(2**curr_norm_abs.loc[most_cut].mean()-1)*100
                                 ,(2**curr_norm_abs.loc[least_cut].mean()-1)*100,pd.Series(curr_tf_col),curr_timepoints],axis=1)
            curr_df.columns = wanted_cols # setting column names to fit the combined final data frame
            fin_df = pd.concat([fin_df,curr_df])  # concatanating curr experiment df and final df

    if len(fin_df) != 0:
        fin_df['Most-Least'] = fin_df['Most']-fin_df['Least'] #new most minus least column
        sns.set(font_scale=1)
        fig,axs = plt.subplots(1,2,figsize=(16,8)) #open figure
        plt.suptitle(lib_info['gene']+':'+str(lib_num))
        # left plot
        curr_ax = 0
        # scattering results
        sns.scatterplot(data=fin_df,x='Most',y='Least',hue='tf',size='timepoint',sizes=(50,200),ax=axs.ravel()[curr_ax],palette='Paired')
        # adding x=0 y=0 and x=y lines
        axs.ravel()[curr_ax].legend([],[], frameon=False)  # removing legend from left plot
        axs.ravel()[curr_ax].axhline(0,c='k')
        axs.ravel()[curr_ax].axvline(0,c='k')
        axs.ravel()[curr_ax].axline((0, 0), slope=1,c='k')
        # ax.set_ylabel('Least cut seqs change in %', fontsize=18)
        axs.ravel()[curr_ax].set_xlabel('Most cut seqs change %', fontsize=18)
        axs.ravel()[curr_ax].set_ylabel('Least cut seqs change %', fontsize=18)
        # right plot (y axis showing most-least)
        curr_ax = 1
        sns.scatterplot(data=fin_df,x='Most',y='Most-Least',hue='tf',size='timepoint',sizes=(50,200),ax=axs.ravel()[curr_ax],palette='Paired')
        # adding x=0 y=0 and x=y lines
        axs.ravel()[curr_ax].axhline(0,c='k')
        axs.ravel()[curr_ax].axvline(0,c='k')
        axs.ravel()[curr_ax].axline((0, 0), slope=1,c='k')
        axs.ravel()[curr_ax].set_xlabel('Most cut seqs change %', fontsize=18)
        axs.ravel()[curr_ax].set_ylabel('Most-Least %', fontsize=18)
        sns.move_legend(axs.ravel()[curr_ax], "upper left", bbox_to_anchor=(1, 1))
        return fig

if __name__ == '__main__':
    x = 1