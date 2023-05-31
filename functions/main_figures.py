import sys
import os
dir_path = '/home/labs/barkailab/tamarj/TF_combinatorics/'
sys.path.insert(0,dir_path)
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pandas as pd
import numpy as np
from functions import parsing
from functions import prom_analysis
from functions import data_processing as dp
from functions import general_data_presentation as gdp
from functions import general_functions as gf
from functions import lib_table_operations as lto
from functions import params
from functions import general_data_presentation as gdp
from functions import basic_figures as bf
from collections import OrderedDict
import re


def heatmaps_tf_lib(tf, lib_num, exp_num, out_path=None):
            lib_info = dp.get_lib_info(lib_num)
            lib_name = lib_info['gene']
            only_intact_pos_change = lto.only_pos_intact_change_calc(tf, lib_num, exp_num, norm_to=None)
            only_mut_pos_change = lto.only_pos_mut_change_calc(tf, lib_num, exp_num, norm_to=None)
            fig, axes = plt.subplots(nrows=2, ncols=1)
            fig.set_size_inches(4, 6)
            col_map = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True)

            plt.subplots_adjust(hspace=3)
            ax = sns.heatmap(only_intact_pos_change * -100, cmap=col_map, annot=True, fmt=".2f",
                             annot_kws={"size": 18.5 / (np.sqrt(len(only_intact_pos_change)) * 2)}, ax=axes[0])
            ax.set_title('Only pos wt change')
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

            ax = sns.heatmap(only_mut_pos_change * -100, cmap=col_map, annot=True, fmt=".2f",
                             annot_kws={"size": 18.5 / (np.sqrt(len(only_mut_pos_change)) * 2)}, ax=axes[1])
            ax.set_title('Only pos mutated change')
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

            plt.suptitle(tf + ' - ' + lib_name)
            if out_path != None:
                plt.savefig(os.path.join(out_path, tf + '_' + lib_name + '_' + 'heatmap.pdf'), bbox_inches='tight')


def seq_scatter_by_tf_motifs(tf, lib_num, exp_num, ax=None, out_path=None, positions=None):
    matplotlib.rcParams['pdf.fonttype'] = 42
    all_prom_df = pd.DataFrame()
    lib_info = dp.get_lib_info(lib_num)
    lib_name = lib_info['gene']
    file_name = tf + '_' + str(lib_num) + '_' + str(exp_num) + '.csv'  # find results file ('tf_libnum_exp.csv')
    res_table = (pd.read_csv(os.path.join(params.RES_PATH, file_name), index_col=0))
    norm_res = dp.norm_reads(res_table)
    sample_filt_norm, _ = dp.rm_samples(tf, lib_num, exp_num, norm_res)
    sample_info = dp.get_samp_info(sample_filt_norm)
    log2_norm = dp.res_log2(sample_filt_norm)
    norm_tp_0 = dp.norm_to_tp_0(log2_norm, sample_info)
    norm_uncut = lto.norm_non_cut(lib_num, norm_tp_0, avg_over=None)
    res_table = lto.convert_fc_to_occ(norm_uncut) * 100

    sem_repeats = dp.sem_bio_rep_tps(res_table, sample_info).iloc[:, 1:]
    mean_repeats = dp.mean_bio_rep_tps(res_table, sample_info).iloc[:, 1:]
    cols = mean_repeats.columns

    mean_repeats['mot_count'] = lto.seq_mut_wt_cmp(mean_repeats, lib_info, tf)[1]
    mean_repeats['Promoter'] = [lib_name] * len(mean_repeats)
    all_prom_df = all_prom_df.append(mean_repeats)
    col_map = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=False, n_colors=len(all_prom_df['mot_count'].unique()))

    if ax ==None:
        fig, ax = plt.subplots(figsize=(6, 4))
    sns.scatterplot(data=all_prom_df, x=cols[0], y=cols[1], hue='mot_count', palette=col_map, edgecolor='none',
                    s=80, zorder=20, ax=ax)
    ax.errorbar(mean_repeats[cols[0]], mean_repeats[cols[1]], xerr=sem_repeats[cols[0]],
                yerr=sem_repeats[cols[1]], ecolor='grey', elinewidth=0.5, color='none', fmt='o')

    wt_seq = ''.join(dp.get_lib_info(lib_num)['wt_at_var_loc'])
    mut_seq = lto.get_mut_var(lib_num)
    extremes = pd.concat([all_prom_df.loc[wt_seq].to_frame().T, all_prom_df.loc[mut_seq].to_frame().T])
    extremes['colors'] = np.array([0,1])
    sns.scatterplot(data=extremes, x=cols[0], y=cols[1], color='red', edgecolor='none',zorder=20, s=80, ax=ax)
    # print(extremes)

    if positions != None:
        tf_pos = lto.get_tf_mot_pos(tf, lib_num)
        wt_var_seq = ''.join(dp.get_lib_info(lib_num)['wt_at_var_loc'])
        non_tf_pos = np.delete(np.array(range(0, len(wt_var_seq))), tf_pos)
        only_pos_wt = lto.get_only_pos_intact(tf_pos, wt_var_seq,mean_repeats)
        only_pos_mut = lto.get_only_pos_intact(non_tf_pos, wt_var_seq, mean_repeats)
        pos_df = all_prom_df.iloc[only_pos_wt,:].to_frame().T
        # print(pos_df)
        sns.scatterplot(data=pos_df, x=cols[0], y=cols[1], color='blue', edgecolor='none', s=80, zorder=20, ax=ax)
        pos_df = all_prom_df.iloc[only_pos_mut,:].to_frame().T
        sns.scatterplot(data=pos_df, x=cols[0], y=cols[1], color='orange', edgecolor='none',  s=80,  zorder=20, ax=ax)

    ax.set_aspect('equal', adjustable='box')
    ax.set_title(tf + ' - '+ lib_name)
    plt.legend(title = tf + ' intact motifs')
    ax.set_xlabel(cols[0])
    ax.set_ylabel(cols[1])
    min_lim = min(ax.get_xlim()[0], ax.get_ylim()[0])
    max_lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.set_xlim(min_lim,max_lim)
    ax.set_ylim(min_lim,max_lim)

    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1.1))
    plt.show()
    if out_path != None:
        fig.savefig(os.path.join(out_path, tf + '_' + lib_name + '_' + 'mot_num_scatter.pdf'), bbox_inches='tight')
    # plt.close()



def heatmaps_by_tf(tf, lib_num, exp_num, tps, out_path=None):
        lib_info = dp.get_lib_info(lib_num)
        lib_name = lib_info['gene']
        only_intact_pos_change, _ = lto.only_pos_intact_change_calc(tf, lib_num, exp_num, norm_to=None, occ=True)
        only_intact_pos_change = only_intact_pos_change.loc[tps]
        only_mut_pos_change = lto.only_pos_mut_change_calc(tf, lib_num, exp_num, None, occ=True).loc[tps]
        fig, axes = plt.subplots(nrows=2, ncols=1)
        fig.set_size_inches(4, 6)
        col_map = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True)

        plt.subplots_adjust(hspace=3)
        ax = sns.heatmap(only_intact_pos_change*100, cmap=col_map, annot=True, fmt=".2f",
                         annot_kws={"size": 18.5 / (np.sqrt(len(only_intact_pos_change)) * 2)}, ax=axes[0])
        ax.set_title('Only pos wt change')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

        ax = sns.heatmap(only_mut_pos_change*100, cmap=col_map, annot=True, fmt=".2f",
                         annot_kws={"size": 18.5 / (np.sqrt(len(only_mut_pos_change)) * 2)}, ax=axes[1])
        ax.set_title('Only pos mutated change')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

        plt.suptitle(tf + ' - ' + lib_name)

        if out_path != None:
            plt.savefig(os.path.join(out_path, tf + '_' + lib_name + '_' + 'heatmap.pdf'), bbox_inches='tight')
        # plt.close()








def obs_exp_calc(tf, lib_num, exp_num, relevant_pos, relevant_combs, tp):
    norm_df = dp.mean_over_bio_reps(lto.norm_res_data(tf, lib_num, exp_num, norm_to=None))
    seqs = list(norm_df.index)
    all_seqs_split = np.array([list(seq) for seq in seqs])
    lib_info = dp.get_lib_info(lib_num)  # get library info
    wt_var_seq = ''.join(lib_info['wt_at_var_loc'])  # get wild-type sequence
    mut_var_seq = lto.get_mut_var(lib_num)
    only_intact_pos_change = lto.only_pos_intact_change_calc(tf, lib_num, exp_num, norm_to=None)
    only_intact_pos_change = only_intact_pos_change.iloc[only_intact_pos_change.index == tp, :]
    expected_dict = {}
    observed_dict = {}

    # calculate observed values:
    for pos in relevant_pos:
        observed_dict[str(pos)] = only_intact_pos_change.values[0][pos]

    for comb in relevant_combs:
        non_comb = np.delete(np.array(range(0, len(wt_var_seq))), comb)  # all positions other than relevant positions
        wt_comb = np.array(list(wt_var_seq))[comb]  # wt seq at comb positions
        mut_non_comb = np.array(list(mut_var_seq))[non_comb]  # mut seq at non-comb positions
        id_rel_wt = np.sum(all_seqs_split[:, comb] == wt_comb, axis=1) == len(comb)
        id_rel_mut = np.sum(all_seqs_split[:, non_comb] == mut_non_comb, axis=1) == len(
            non_comb)  # ids of seqs with mut seq at non-comb positions
        observed_dict[''.join(map(str, comb))] = norm_df.loc[np.array(seqs)[id_rel_wt & id_rel_mut], :][tp].values[0]

    # calculate observed expected:
    for comb in relevant_combs:
        sub_comb = gf.get_pair_combs(comb, len(comb) - 1)
        sub_comb_observed = []
        for curr_sub_comb in sub_comb:
            curr_key = ''.join(map(str, curr_sub_comb))
            sub_comb_observed.append(observed_dict[curr_key])
        expected_dict[''.join(map(str, comb))] = (1 - np.prod(1 - (-1 * (np.array(sub_comb_observed))))) * 100

    obs_comb_ids = np.array([i for i, key in enumerate(observed_dict.keys()) if len(key) > 1])
    # print(expected_dict, observed_dict)
    expected = np.array(list((expected_dict.values())))
    observed = np.array(list((observed_dict.values())))[obs_comb_ids] * -100

    return expected, observed


def plot_obs_exp(tf_list, tp, threshold=None, relevant_tfs=None, tf_mot_num=None, out_path=None):
    # matplotlib.rcParams['pdf.fonttype'] = 42
    col_names = ['Observed', 'Expected', 'Obs-Exp', 'Motif_combination', 'Relevant TF motifs', 'TF', 'Library',
                 'Library_number']
    df_obs_exp = pd.DataFrame(columns=col_names)
    for tf in tf_list:
        tf_files = dp.find_tf_files(tf)
        for file in tf_files:
            _, lib_num, exp_num = file.split('.')[0].split('_')
            if exp_num == '26':
                continue
            lib_info = dp.get_lib_info(lib_num)  # get library info
            lib_name = lib_info['gene']
            tf_mot_num_in_lib = len(lto.get_tf_mot_pos(tf, lib_num))

            if ((tf_mot_num != None) and (tf_mot_num != tf_mot_num_in_lib)):
                continue

            only_intact_pos_change = lto.only_pos_intact_change_calc(tf, lib_num, exp_num, norm_to=None)

            tp_cond = only_intact_pos_change.index == tp
            only_intact_pos_change = only_intact_pos_change.iloc[tp_cond, :]

            if threshold != None:
                cond_threshold = only_intact_pos_change < threshold
                cond_threshold = cond_threshold.values[0]
            else:
                cond_threshold = ~np.ones(np.shape(only_intact_pos_change), dtype=bool)[0]

            if relevant_tfs != None:
                relevant_tfs_pos = []
                for curr_tf in relevant_tfs:
                    tf_name_split = re.split(r'(\d+)', curr_tf)[0]
                    tf_pos = [i for i, mot in enumerate(lib_info['mut_by_tf']) if tf_name_split in mot]
                    relevant_tfs_pos = relevant_tfs_pos + tf_pos
                cols_ids = cond_threshold
                cols_ids[np.array(relevant_tfs_pos)] = True
            else:
                cols_ids = cond_threshold

            relevant_pos = np.where(cols_ids)[0]
            if len(relevant_pos) == 0:
                continue

            relevant_combs = []
            perm_sizes = list(range(2, len(relevant_pos) + 1))  # sizes of possible combinations
            for curr_perm_size in perm_sizes:  # iterate over combination sizes
                relevant_combs = relevant_combs + gf.get_pair_combs(relevant_pos, curr_perm_size)

            if len(relevant_combs) == 0:
                continue

            expected, observed = obs_exp_calc(tf, lib_num, exp_num, relevant_pos, relevant_combs, tp)
            relevant_tfs_by_pos = only_intact_pos_change.columns[relevant_pos].values
            curr_df = (pd.DataFrame([observed, expected, observed - expected, relevant_combs,
                                     np.array([[relevant_tfs_by_pos], ] * len(expected)),
                                     np.repeat(tf, len(expected)), np.repeat(lib_name, len(expected)),
                                     np.repeat(int(lib_num), len(expected))]).T)
            curr_df.columns = col_names
            df_obs_exp = df_obs_exp.append(curr_df)

    tf_mot_num = []
    only_one_tf = []
    for i in range(len(df_obs_exp['Relevant TF motifs'])):
        curr_comb_pos = df_obs_exp['Motif_combination'].values[i]
        curr_lib_num = df_obs_exp['Library_number'].values[i]
        curr_lib_info = dp.get_lib_info(curr_lib_num)
        tf_list = np.array(curr_lib_info['mut_by_tf'])[curr_comb_pos]
        tf_name = df_obs_exp['TF'].values[i]
        tf_name_split = re.split(r'(\d+)', tf_name)[0]
        tf_mot_num.append(len([s for s in tf_list if tf_name_split in s]))
        uniqe_tf_comb = np.unique(tf_list)
        if len(uniqe_tf_comb) < 2:
            only_one_tf.append(False)
        else:
            only_one_tf.append(True)
    df_obs_exp['TF_motif_number'] = tf_mot_num
    df_obs_exp['TF_only'] = only_one_tf

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))
    ax1 = sns.scatterplot(data=df_obs_exp, x='Observed', y='Obs-Exp', hue='Library', style='TF', size='TF_motif_number',
                          sizes=(100, 350), ax=axes[0])
    ax1.legend(loc='upper left', bbox_to_anchor=(1, 1.1))
    max_lim = max(ax1.get_xlim()[1], ax1.get_ylim()[1])
    min_lim = min(ax1.get_xlim()[0], ax1.get_ylim()[0])
    # ax1.set_xlim(min_lim, max_lim)
    # ax1.set_ylim(min_lim, max_lim)
    ax1.plot(ax1.get_xlim(), (0, 0), ls="--", c=".3")
    ax1.set_aspect('equal', 'box')
    comb_sizes = np.array([len(row) for row in df_obs_exp['Motif_combination'].values])
    df_obs_exp_comb_2 = df_obs_exp.iloc[comb_sizes == 2, :]
    pair_combs = df_obs_exp_comb_2['Motif_combination'].values
    mot_dist = []
    for i, comb in enumerate(pair_combs):
        curr_mut_loc = np.array(dp.get_lib_info(df_obs_exp_comb_2['Library_number'].values[i])['mut_loc'])
        mot_dist.append(np.diff(curr_mut_loc[comb])[0])
    df_obs_exp_comb_2['Motif distances'] = mot_dist

    obs_thershold = 20
    obs_cond = df_obs_exp_comb_2['Observed'].values > obs_thershold
    df_obs_exp_comb_2_trsh = df_obs_exp_comb_2.iloc[obs_cond, :]

    labels = list(df_obs_exp_comb_2_trsh['Library'].unique())

    sns.scatterplot(data=df_obs_exp_comb_2_trsh, x='Observed', y='Obs-Exp', hue='Library', style='TF',
                          s=200, style_order=list(df_obs_exp_comb_2_trsh['TF'].unique()), hue_order=labels,
                          legend=False, ax=axes[1])

    df_obs_exp_comb_2_one_tf = df_obs_exp_comb_2_trsh.iloc[df_obs_exp_comb_2_trsh['TF_only'].values, :]
    ax2 = sns.scatterplot(data=df_obs_exp_comb_2_one_tf, x='Observed', y='Obs-Exp', hue='Library', style='TF',
                          s=200, edgecolor='black', linewidth=1.2,
                          style_order=list(df_obs_exp_comb_2_trsh['TF'].unique()),
                          hue_order=labels, ax=axes[1])

    ax2.legend(loc='upper left', bbox_to_anchor=(1, 1.1))
    # max_lim = max(ax2.get_xlim()[1], ax.get_ylim()[1])
    # min_lim = min(ax2.get_xlim()[0], ax.get_ylim()[0])
    ylim_vals = ax2.get_ylim()
    xlim_vals = ax2.get_xlim()

    # ax2.set_xlim(obs_thershold, max_lim)
    # print(-abs(max(ylim_vals)), abs(max(ylim_vals)))
    ax2.set_ylim(-abs(max(ylim_vals)), abs(max(ylim_vals)))
    ax2.plot(ax2.get_xlim(), (0, 0), ls="--", c=".3")
    ax2.set_xlim(xlim_vals[0], xlim_vals[1])
    # plt.axis('square')

    if out_path != None:
        plt.savefig(os.path.join(out_path, 'obs_exp_threshold.pdf'), bbox_inches='tight')

    return df_obs_exp_comb_2_one_tf

