import numpy as np
import scipy.io as sio
from functions import parsing
from Bio.Seq import Seq


def get_promoter_pos(gene_ind, promoter_length=None, into_gene=None):
    ''' This function gets gene id (with optional parameters of promoter length and number of bs into the gene body)
    and returns the genomic location (pos 0 = chr number, pos 1 = start pos 2 = end start pos 3 = strand).'''
    pos = STEIN[gene_ind, :]
    chr_num = pos[0]
    strand = np.sign(pos[2]-pos[1])
    if promoter_length is None:
        promoter_length_vec = PROMOTER_LENGTH
    else:
        promoter_length_vec = np.zeros(shape=np.shape(STEIN)[0])
        promoter_length_vec[~np.isnan(STEIN[0:, 0])] = promoter_length
    curr_promoter_length = promoter_length_vec[gene_ind]
    if strand > 0:
        start_pos = pos[1] - curr_promoter_length
        end_pos = pos[1]
        if into_gene != None:
            end_pos += into_gene
    else:
        start_pos = pos[1] + 1
        end_pos = start_pos + curr_promoter_length
        if into_gene != None:
            start_pos -= into_gene
    #assert(end_pos - start_pos == curr_promoter_length)
    return np.array([chr_num, start_pos, end_pos, strand]).astype(int)



def get_seq_by_pos(pos):
    '''This function gets genomic positions (pos 0 = chr number, pos 1 = start pos 2 = end start pos 3 = strand) and returnt
    the sequennce and the reverse complement sequence.'''
    chr_num = "chr" + str(int(pos[0]+1))
    if pos[3] > 0:
        seq = Cer_Genome[chr_num][pos[1]:pos[2]]
        rc_seq = str(Seq(seq).reverse_complement())
    else:
        rc_seq = Cer_Genome[chr_num][pos[1]:pos[2]]
        seq = str(Seq(rc_seq).reverse_complement())
    return [seq, rc_seq]


if __name__ == '__main__':
    x = 1