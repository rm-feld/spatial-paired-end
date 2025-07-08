import numpy as np
from io import StringIO
from Bio import Phylo
import subprocess 

# ---------------------------------------------
# high level path info/utils
REF = "/Users/rfeld/Documents/Research/SPATIAL/spatial_24/ref/chr21.fa"
REFNAME = "chr21"
OUTDIR = "/Users/rfeld/Documents/Research/SPATIAL/output"
# READ_PREFIX = "h"
# SEQ_PREFIX = "h"

BUFFER = 20
# ---------------------------------------------


REFNAME = "/Users/rfeld/Documents/Research/SPATIAL/spatial_24/ref/chr21.fa"

def mutate(b_rel_idx, b_rel_start, b_length, copies):
    region = copies * b_rel_idx[b_rel_start:b_rel_start + b_length]
    return b_rel_idx[:b_rel_start] + region + b_rel_idx[b_rel_start + b_length:]

def generate_random_mutation(b_rel_idx, gdir):
    num_regions = len(b_rel_idx)
    CN_size = max(min(round(np.random.exponential(gdir["cn_length_mean"] / gdir["resolution"])), num_regions), 1)
    CN_type = np.random.binomial(1, 0.5)
    CN_copies = np.random.geometric(0.5) + 1 if CN_type == 1 else 0
    b_start_idx = np.random.randint(num_regions - CN_size + 1)

    return {
        "copies": CN_copies,
        "start_idx": b_start_idx,
        "length_blocks": CN_size,
        "region_length": CN_size * gdir["resolution"],
        "total_length": CN_copies * CN_size * gdir["resolution"],
        "prev_start_idx": b_start_idx * gdir["resolution"] + 1,
        "prev_end_idx": (b_start_idx + CN_size) * gdir["resolution"]
    }

def apply_mutation_and_trace(blocks, mutation):
    return mutate(blocks, mutation["start_idx"], mutation["length_blocks"], mutation["copies"])

def simulate_cell_mutations(blocks, n_mut, gdir):
    mutations = []
    for _ in range(n_mut):
        md = generate_random_mutation(blocks, gdir)
        blocks = apply_mutation_and_trace(blocks, md)
        mutations.append(md)
    return blocks, mutations

def log_mutation(logf, name, allele, mut_id, mutation, gdir):
    # Log mutation details in a TSV format
    logf.write(f"{name}\t{allele}\t{mut_id}\t"
               f"{mutation['start_idx']}\t"
               f"{mutation['start_idx'] + mutation['length_blocks'] - 1}\t"
               f"{mutation['prev_start_idx']}\t"
               f"{mutation['prev_end_idx']}\t"
               f"{mutation['copies']}\t"
               f"{mutation['length_blocks']}\t"
               f"{mutation['region_length']}\n")

def simulate_diploid_cell(parent_blocks, n0, n1, inherited_mutations, logf, name, gdir):
    all_mutations = {}
    final_blocks = {}
    for allele, n_mut in [(0, n0), (1, n1)]:
        blocks = parent_blocks[allele]
        inherited = inherited_mutations[allele].copy()
        
        logf.write(f"\n{name} allele {allele}: Starting with {len(blocks)} blocks\n")
        
        # Apply inherited mutations
        for md in inherited:
            blocks = apply_mutation_and_trace(blocks, md)
        
        # Generate new mutations
        blocks, new_muts = simulate_cell_mutations(blocks, n_mut, gdir)
        all_mutations[allele] = inherited + new_muts
        final_blocks[allele] = blocks

        logf.write(f"{name} allele {allele}: {n_mut} new mutations, final length: {len(blocks)}\n")
        
        # Log all mutations with detailed info
        for i, md in enumerate(new_muts):
            log_mutation(logf, name, allele, i + 1, md, gdir)

    return all_mutations, final_blocks

def logme(log):
    if log == "console":
        return print
    elif log:
        def _logme(s):
            with open(log, "a+") as f:
                f.write(s) 
        return _logme
    else: 
        return lambda x: None

def construct(nwk_string, diploid_blocks, gdir, log="console"):
    # Setup logger
    logger = logme(log)

    # Open the log file if needed
    with open(log, "w") if log != "console" else None as logf:
        logf.write("Name\tAllele\tMutation_ID\tStart_Block\tEnd_Block\tStart_BP\tEnd_BP\tCopies\tLength_Blocks\tLength_BP\n")

        tree = Phylo.read(StringIO(nwk_string), "newick")
        tree.ladderize()
        parent_map = {child: clade for clade in tree.find_clades(order="level") for child in clade.clades}
        cumulative_genomes = {}
        mutation_records = {}

        for node in tree.find_clades(order="level"):
            if not node.name:
                continue

            logger(f"\nProcessing {node.name}...")
            parent = parent_map.get(node)

            if parent is None:
                cumulative_genomes[node.name] = {0: diploid_blocks.copy(), 1: diploid_blocks.copy()}
                mutation_records[node.name] = {0: [], 1: []}
                continue
            
            inherited_mutations = mutation_records[parent.name]
            parent_blocks = {allele: diploid_blocks.copy() for allele in [0, 1]}

            if gdir["num_mutations"] == "random":
                n0 = np.random.choice(gdir["allowed_num"], p=gdir["num_mutations_w"])
                n1 = np.random.choice(gdir["allowed_num"], p=gdir["num_mutations_w"])
            else:
                n0 = gdir["num_mutations"][node.name][0]
                n1 = gdir["num_mutations"][node.name][1]
                
            muts, final_blocks = simulate_diploid_cell(
                parent_blocks, n0, n1, inherited_mutations, logf, node.name, gdir
            )

            mutation_records[node.name] = muts
            cumulative_genomes[node.name] = final_blocks

    return tree, mutation_records, cumulative_genomes

def write_fasta(block_index, header, ref, fname, fpath):
    seq = [ref[i] for i in block_index]
    fasta_loc = f"{fpath}/{fname}.fa"
    with open(fasta_loc, "w") as f:
        f.write(f">{header}\n")
        for block in seq:
            for i in range(0, len(block), 60):
                f.write(block[i:i+60] + "\n")
    
    return fasta_loc

def write_fasta_slice(block_index, ref, left, right, buffer, 
                      header, name, fpath, gdir):
    # ASSUME BLOCK INDICES FOR left, right
    l = max(0, left - buffer)
    r = min(right + buffer, len(block_index))

    left = l * gdir["resolution"] + 1
    right = r * gdir["resolution"]

    fasta_loc = f"{fpath}/{name}_{left}_{right}.fa"

    seq = [ref[i] for i in block_index][l:r]
    with open(fasta_loc, "w") as f:
        f.write(f">{header}\n")
        for block in seq:
            for i in range(0, len(block), 60):
                f.write(block[i:i+60] + "\n")
    
    return fasta_loc

def clean_ref(refname):
    ref = []
    with open(refname, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                ref.append(line.strip())

    # Each line ~60bp → group 16 lines into 960bp block
    ref_blocks = []
    for i in range(0, len(ref), 16):
        block = "".join(ref[i:i+16])
        ref_blocks.append(block)

    # Optional: track index
    ref_idx = list(range(len(ref_blocks)))
    reflength = 960 * (len(ref_blocks) - 1) + len(ref_blocks[-1])

    return ref, ref_idx, reflength

def generate(gdir, outdir, seqdir, readdir, refname = REFNAME, log = "console"):

    cna_call = f"cnasim -m 0 -n {gdir['numcells']} -r1 {REF} -L 45090682 -o {seqdir} -j {gdir['founder_mult']} -k {gdir['resolution']}"
    subprocess.run(cna_call.split())

    with open(f'{seqdir}/tree.nwk', "r") as f:
        nwk_string = f.readline()
    
    ref, ref_idx, reflength = clean_ref(refname) 

    tree, mutation_records, cumulative_genomes = construct(nwk_string, ref, 
                                                           gdir, log=log)
    
    cellnames = list(mutation_records.keys())[1:]
    _idx = np.random.choice(len(cellnames))

    cellname = cellnames[_idx]

    


    
