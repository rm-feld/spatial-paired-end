import numpy as np
from .consts import *

def insertion(seq, m, pos1, pos2):
    start = np.where(seq == pos1)[0][0]
    stop = start + pos2 - pos1 + 1
    region = np.tile(seq[start:stop], m)
    return np.concatenate((seq[:start], region, seq[stop:]), axis = None)

def deletion(seq, pos1, pos2):
    start = np.where(seq == pos1)[0][0]
    stop = start + pos2 - pos1 + 1
    return np.concatenate((seq[:start], seq[stop:]), axis = None)

def random_mutations(seq):
    n = np.random.choice(POP_MUT, p = POP_MUT_WEIGHT)
    smax = READ_SEP[1]
    lmin, lmax = MUT_LENGTH
    mutations = []

    for i in range(n):
        mut_len = np.random.randint(lmin, lmax)
        start_loc = np.random.randint(smax, len(seq) - smax - mut_len)
        end_loc = start_loc + mut_len
        if np.random.randint(2) == 0:
            m = np.random.choice(DUP_NUM, p = DUP_WEIGHT)
            seq = insertion(seq, m, start_loc, end_loc)
            mutations.append((start_loc, end_loc, m))
        else:
            seq = deletion(seq, start_loc, end_loc)
            mutations.append((start_loc, end_loc, 0))
    
    return seq, mutations