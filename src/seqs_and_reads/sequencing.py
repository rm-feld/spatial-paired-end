import numpy as np
from .consts import * 

def paired_end_reads(seq, prop):
    """ generates paired end reads of sequence with intensity moderated by frequency
    of population

    Parameters
    ----------
    seq : Array
        Array of integers representing position to which a read should be mapped
        on reference genome 
    prop : float
        proportion of total population which this genome makes up
    Returns
    -------
    coords : Array
        n x 2 array of points
    """    
    R0, R1 = READ_SEP
    # identify approximate number of reads we need to get expected coverage
    n = int((len(seq) * prop) // (2 * READ_LEN) * COVERAGE)
    l = len(seq)

    def choose_r1s(sep):
        return np.random.randint(0, l - sep - 2 * READ_LEN - 1)
    choose_r1s_v = np.vectorize(choose_r1s)

    # choose separation randomly
    seps = np.random.randint(R0, R1, size = n)
    r1s = choose_r1s_v(seps) # pick first read loc with above restriction
    r1e = r1s + READ_LEN
    r2s = r1e + seps
    r2e = r2s + READ_LEN

    # return coords as an n x 2 array of coordindate values
    coords = np.asarray([seq[r1s], seq[r2e]]).T
    
    return coords

def find_coverage(reads, ref):
    #TODO: update to take only coords
    l = len(ref)
    reads = np.concatenate([np.concatenate(read).flat for read in reads])
    return np.bincount(reads, minlength = l)