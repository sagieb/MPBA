import numpy as np
import scipy.io as sio
from functions import parsing
from Bio.Seq import Seq

STEIN = np.load('/home/labs/barkailab/tamarj/TF_combinatorics/data/genome_annot/stein.npy')
Cer_Genome = parsing.decompress_pickle('/home/labs/barkailab/tamarj/TF_combinatorics/data/genome_annot/cerGenome/cerGenome.pbz2')
Promoter_Length = np.load('/home/labs/barkailab/tamarj/TF_combinatorics/data/genome_annot/prompter_length_GP.npy')
GENE_NAMES = sio.loadmat('/home/labs/barkailab/tamarj/TF_combinatorics/data/genome_annot/gene_names.mat', simplify_cells=True)['gene_names']

#This function gets gene id (with optional parameters of promoter length and number of bs into the gene body)
# and returns the genomic location (pos 0 = chr number, pos 1 = start pos 2 = end start pos 3 = strand)
def get_promoter_pos(gene_ind, promoter_length=None, into_gene=None):
    pos = STEIN[gene_ind, :]
    chr_num = pos[0]
    strand = np.sign(pos[2]-pos[1])
    if promoter_length is None:
        promoter_length_vec = Promoter_Length
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


#This function gets genomic positions (pos 0 = chr number, pos 1 = start pos 2 = end start pos 3 = strand) and returnt
# the sequennce and the reverse complement sequence
def get_seq_by_pos(pos):
    chr_num = "chr" + str(int(pos[0]+1))
    if pos[3] > 0:
        seq = Cer_Genome[chr_num][pos[1]:pos[2]]
        rc_seq = str(Seq(seq).reverse_complement())
    else:
        rc_seq = Cer_Genome[chr_num][pos[1]:pos[2]]
        seq = str(Seq(rc_seq).reverse_complement())
    return [seq, rc_seq]
