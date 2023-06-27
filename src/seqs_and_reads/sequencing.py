import numpy as np
from .consts import * 

def paired_end_reads(seq, prop):
    R0, R1 = READ_SEP
    # identify approximate number of reads we need to get expected coverage
    n = (len(seq) * prop) // (2 * READ_LEN) * COVERAGE
    l = len(seq)

    # initialize reads and ends
    reads = []
    coords = []

    print("expected iterations:", n)
    for i in range(n):
        # choose separation randomly 
        sep = np.random.randint(R0, R1)
        # pick first read loc with above restriction
        r1s = np.random.randint(0, l - sep - 2 * READ_LEN - 1)
        r1e = r1s + READ_LEN
        r2s = r1e + sep 
        r2e = r2s + READ_LEN
        reads.append([seq[r1s:r1e], seq[r2s:r2e]])
        coords.append((seq[r1s], seq[r2e]))
        if i % 40 == 0:
            print(f"{n - i} iterations to go!")
    
    return reads, coords

def find_coverage(reads, ref):
    l = len(ref)
    reads = np.concatenate([np.concatenate(read).flat for read in reads])
    return np.bincount(reads, minlength = l)