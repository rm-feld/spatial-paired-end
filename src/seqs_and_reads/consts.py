

# GENOME SPECIFICATIONS
FRAGMENT_LENGTH = [50000, 1000000] # size of fragment being sequenced

# MUTATION SPECIFICATIONS
MAX_POPS = 7
POP_MUT = [1, 2, 3] # number of distinct mutations added per population
POP_MUT_WEIGHT = [0.2, 0.5, 0.3] # population frequencies
MUT_LENGTH = [4000, 15000]
DUP_NUM = [2, 3, 4]
DUP_WEIGHT = [0.6, 0.3, 0.1]

# READ SPECIFICATIONS
TOL = 20 # number of 
READ_SEP = [1000, 3000] # range of separations of reads
READ_LEN = 150
COVERAGE = 30

# EVALUATION SPECS
ETOL = 0.05 # frequency diff at which populations are considered 'the same'