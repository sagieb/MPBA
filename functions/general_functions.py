import itertools
from itertools import combinations
from itertools import permutations
import seaborn as sns
import numpy as np

# This function generates all the unique pairwise combinations of elements in a vectors:
def get_pair_combs(vec, comb_size):
    all_combs = [list((map(int, comb))) for comb in combinations(vec, comb_size)]
    return all_combs


def get_substring(my_string, wanted_ids):
    # this function gets a string and returns a substring based on wanted ids
    # for example get_substring('Sagikie',[0,1,2,3,6]) gives you back the substring 'Sagie'
    substring = ''.join([sub_str for i, sub_str in enumerate(my_string) if i in wanted_ids])
    return substring


def list_substring(my_list,wanted_ids):
    # this functionr reciveves a list of strings and returns a list containing the substring indicated by wanted_ids for each string in the list
    # this function uses the "get_substring" function
    sub_seqs = [get_substring(seq,wanted_ids) for seq in my_list]
    return sub_seqs


def plot_arc(three_points_cords, ax, color):
    # Linear length along the line:
    three_points_cords = np.array(three_points_cords)
    distance = np.cumsum( np.sqrt(np.sum( np.diff(three_points_cords, axis=0)**2, axis=1 )) )
    distance = np.insert(distance, 0, 0)/distance[-1]
    # Build a list of the spline function, one for each dimension:
    splines = [UnivariateSpline(distance, coords, k=2, s=.2) for coords in three_points_cords.T]
    # Computed the spline for the asked distances:
    alpha = np.linspace(0, 1, 75)
    points_fitted = np.vstack(spl(alpha) for spl in splines).T
    # Graph:
    ax.plot(*three_points_cords.T, 'ok')
    ax.plot(*points_fitted.T, c=color)


def generate_color_dict(palette_name, keys, jump):
    '''This function gets a name of a seaborn color palette and an array with key names and generates
    a dictionary of colors.'''
    colors = list(sns.color_palette(palette_name))
    colors = [colors[i] for i in np.arange(0, len(keys) * jump, jump)]
    colors = [np.array(list(c[0:3])) for c in colors]
    color_dict = dict(zip(keys, colors))
    return color_dict


def comb_2_arrays(list_1, list_2):
    permut = itertools.permutations(list_1, len(list_2))
    all_perms = []
    for comb in permut:
        all_perms += (list(zip(comb, list_2)))
    return all_perms


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

if __name__ == '__main__':
    x = 1

