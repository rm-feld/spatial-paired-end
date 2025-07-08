import subprocess
import numpy as np 
import matplotlib.pyplot as plt 
import csv 
import tree as t
from functools import reduce
import pandas as pd
import copy
import argparse
import random
import gzip
import os


# ---------------------------------------------
# high level path info/utils
REF = "/Users/rfeld/Documents/Research/SPATIAL/spatial_24/ref/chr21.fa"
REFNAME = "chr21"
OUTDIR = "/Users/rfeld/Documents/Research/SPATIAL/output"
# READ_PREFIX = "h"
# SEQ_PREFIX = "h"

BUFFER = 100000 #CHECK: should be sufficient given size of reads

samtools = "/Users/rfeld/Documents/Research/SPATIAL/spatial_24/DWGSIM/samtools/samtools"
# ---------------------------------------------

def _save_settings(expdir, seq_dict, read_dict, exp_dict):
    with open(f"{expdir}/log.txt", 'w+') as f:
        f.write("seq_dict = {\n")
        for key, value in seq_dict.items():
            f.write(f"\t{key}: {value},\n")
        f.write("}\n")

        f.write("\nread_dict = {\n")
        for key, value in read_dict.items():
            f.write(f"\t{key}: {value},\n")
        f.write("}\n")

        f.write("\nexp_dict = {\n")
        for key, value in exp_dict.items():
            f.write(f"\t{key}: {value},\n")
        f.write("}\n\n")

    print(f"Experiment log at {expdir}/log.txt")


def _parse_ref(d):
    with open(REF, "r") as f:
        ref = []
        for row in f:
            ref.append(row.strip())
    ref = ref[1:]

    # join rows by resolution
    _rows_per_res = d["resolution"] // len(ref[0])
    _num_groups = len(ref) // _rows_per_res 
    _left_groups = len(ref) % _rows_per_res

    _ref = []
    for i in range(_num_groups):
        _ref.append(''.join(ref[i * _rows_per_res : (i + 1) * _rows_per_res]))
    if _left_groups != 0:
        _ref.append(''.join(ref[-_left_groups:]))

    seqlen = len(_ref) * d["resolution"]

    return seqlen, _ref

def _writefa(filedir, seq, cellname, allele, d, targets = []):
    with open(f'{filedir}/{cellname}.{allele}.fa', "w") as f:
        f.write(f">{cellname}.{allele}\n")
        for line in seq:
            f.write(line + "\n")
    
    if len(targets) > 0:
        for target in targets:
            left = target[0] // d["resolution"] ; right = target[1] // d["resolution"]
            l = left * d["resolution"]; r = right * d["resolution"] 

            augmented_filepath = f'{filedir}/{cellname}.{allele}.{l}.{r}.fa'
            with open(augmented_filepath, "w") as f:
                f.write(f">{cellname}.{allele}.{l}.{r}\n")
                for line in seq[left:right]:
                    f.write(line + "\n")
    return f'{filedir}/{cellname}.{allele}.fa', augmented_filepath

def _recursive_traversal(t):
    if len(t) == 0:
        return []
    elif type(t[-1]) == str:
        return [t[-1]] + _recursive_traversal(t[:-1])
    else:
        return _recursive_traversal(t[-1]) + _recursive_traversal(t[:-1])
    
def _read_focal(fname, cells):
    # cell names 
    out = dict([(c, []) for c in cells])
    out["founder"] = []
    # read in
    raw_focal = []
    with open(fname, "r") as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            raw_focal.append(row)

    for row in raw_focal[1:]:
        out[row[0]].append(row[2:])
    return out

def _compute_lr(mode):
    if mode == 0:
        pass
    elif mode == 1:
        def compute_lr(row):
            if row["strand"] == 1:
                l = row["pos"]; r = row["mpos"]
            else:
                l = row["mpos"]; r = row["pos"]
            return pd.Series({"l": l, 
                              "r": r})
    elif mode == 2:
        pass

    return compute_lr

def _get_read_names(file_path, cellname):
    read_names = set()
    with gzip.open(file_path, 'rt') as f:
        while True:
            header = f.readline().strip()
            f.readline()  # seq
            f.readline()  # plus
            f.readline()  # qual
            if not header:
                break
            if header.startswith(f"@{cellname}"):  # Mapped reads
                read_name = header.split()[0].replace("/1", "").replace("/2", "")
                read_names.add(read_name)
    return read_names

def _filter_reads(input_file, output_file, names_to_keep):
    with gzip.open(input_file, 'rt') as f_in, open(output_file, 'w') as f_out:
        while True:
            header = f_in.readline()
            seq = f_in.readline()
            plus = f_in.readline()
            qual = f_in.readline()
            if not header:
                break
            read_name = header.strip().split()[0].replace("/1", "").replace("/2", "")
            if read_name in names_to_keep:
                f_out.write(header + seq + plus + qual)

def genstep(d, seqdir, only_insertions = False, verbose = False, 
            pick = True, easy = False, no_intersect = True, log = "", seq = True):
    cna_call = f"cnasim -m 0 -n {d['numcells']} -r1 {REF} -L 45090682 -o {seqdir} -j {d['founder_mult']} -k {d['resolution']}"
    subprocess.run(cna_call.split())

    with open(f'{seqdir}/tree.nwk', "r") as f:
        nwk_string = f.readline()

    # flatten for list access 
    traversal = t.parse_newick_structure(nwk_string)
    cells_flattened = _recursive_traversal(traversal)

    focal_loc = f'{seqdir}/focal_events.tsv'
    focal = _read_focal(focal_loc, cells_flattened)

    root, celldict = t.build_tree(traversal, focal, d)

    seqlen, _ref = _parse_ref(d)
    root.ref = _ref

    if pick and seq:
        fp, aug_fp, tm = pickcell(celldict, seqdir, d, easy, no_intersect, 
                                  log = log)
        return fp, aug_fp, tm, root, celldict, focal
    elif pick and not seq:
        tm = pickcell(celldict, seqdir, d, easy, no_intersect, 
                                  log = log, seq = False)
        return tm, root, celldict, focal
    
def populate_fa(tm, celldict, seqdir, d, log = ""):
    tm_cellname = tm.edgeid.split('.')[0]
    tm_allele = tm.allele 
    relloc = (tm.relstart, tm.relstart + tm.length * tm.num_copies)
    #FIXME should probably compute length of mutated sequence
    _left, _right = max(0, relloc[0] - BUFFER), relloc[1] + BUFFER
    tm_seq = celldict[tm_cellname].sequence(tm_allele)
    fp, aug_fp = _writefa(seqdir, tm_seq, tm_cellname, 
                            tm_allele, d, [(_left, _right)])
        
    if log != "":
        with open(log, "a+") as f:
            f.write(f"Chosen mutation specs: seq at {aug_fp}\n")
            f.write(f"\tcell: {tm_cellname}\n")
            f.write(f"\tallele: {tm_allele}\n")
            f.write(f"\tnum_copies: {tm.num_copies}\n")
            f.write(f"\trelstart: {tm.relstart}\n")
            f.write(f"\trefstart: {tm.refstart}\n")
            f.write(f"\tlength: {tm.length}\n\n")
    
    return fp, aug_fp


def pickcell(celldict, seqdir, d, easy = False, no_intersect = True, 
             verbose = False, log = "", seq = True):
    if no_intersect:
        candidates = dict()
        for cell_name in celldict.keys():
            cell = celldict[cell_name]
            no_0 = t.check_non_overlapping(cell, 0)
            no_1 = t.check_non_overlapping(cell, 1)
            non_overlapping = no_1 + no_0 
            candidates[cell_name] = non_overlapping[::]
            flattened_candidates = [m for cellms in candidates.values() for m in cellms]
            sorted_candidates = sorted(flattened_candidates, 
                                       key = lambda m: m.length * m.num_copies)
            insertions = [sc for sc in sorted_candidates if sc.num_copies != 0]
        # For proof of concept, pick smallest region
        if easy:
            tm = insertions[0]
        else:
            _n = random.randint(0, len(sorted_candidates) - 1)
            tm = sorted_candidates[_n]
    else:
        raise NotImplementedError
    
    tm_cellname = tm.edgeid.split('.')[0]
    tm_allele = tm.allele 
    relloc = (tm.relstart, tm.relstart + tm.length * tm.num_copies)
    #FIXME should probably compute length of mutated sequence
    _left, _right = max(0, relloc[0] - BUFFER), relloc[1] + BUFFER
    tm_seq = celldict[tm_cellname].sequence(tm_allele)

    if seq:
        fp, aug_fp = _writefa(seqdir, tm_seq, tm_cellname, 
                            tm_allele, d, [(_left, _right)])
        
        if log != "":
            with open(log, "a+") as f:
                f.write(f"Chosen mutation specs: seq at {aug_fp}\n")
                f.write(f"\tcell: {tm_cellname}\n")
                f.write(f"\tallele: {tm_allele}\n")
                f.write(f"\tnum_copies: {tm.num_copies}\n")
                f.write(f"\trelstart: {tm.relstart}\n")
                f.write(f"\trefstart: {tm.refstart}\n")
                f.write(f"\tlength: {tm.length}\n\n")
    
        return fp, aug_fp, tm
    
    else: 
        return tm

def readstep(d, readdir, rp, tm, aug_fp, verbose = False, log = ""):
    cprint = lambda x: print(x) if verbose else ""
    cellname = tm.edgeid.split(".")[0]

    read1 = f"{readdir}/{rp}.bwa.read1.fastq"
    read2 = f"{readdir}/{rp}.bwa.read2.fastq"

    line =f"dwgsim -N {int(d['numreads'])} -1 {d['readlen']} -2 {d['readlen']} -d {d['outerdist']} -S {d['readmode']} {aug_fp} {rp}"
    subprocess.call(line, shell = True, cwd=readdir)

    compute_lr = _compute_lr(d['readmode'])

    # Define input and output files
    cwd = readdir
    read1 = f"{readdir}/{rp}.bwa.read1.fastq.gz"
    read2 = f"{readdir}/{rp}.bwa.read2.fastq.gz"
    dual_read1 = f"{readdir}/filtered_read1.fastq"
    dual_read2 = f"{readdir}/filtered_read2.fastq"
    unmapped = f"{readdir}/unmapped.fastq"

    # 1. Extract read names from both files
    cprint("Extracting read names...")
    read1_names = _get_read_names(read1, cellname)
    read2_names = _get_read_names(read2, cellname)

    # 2. Identify dually mapped and unmapped reads
    dually_mapped = read1_names & read2_names
    unmapped_reads = (read1_names | read2_names) - dually_mapped

    cprint(f"Found {len(dually_mapped)} dually mapped reads.")
    cprint(f"Found {len(unmapped_reads)} unmapped reads.")

    # 4. Filter dually mapped reads
    cprint("Filtering dually mapped reads...")
    _filter_reads(read1, dual_read1, dually_mapped)
    _filter_reads(read2, dual_read2, dually_mapped)

    # 5. Filter unmapped reads (combine both ends)
    cprint("Filtering unmapped reads...")
    with open(unmapped, 'w') as f_out:
        _filter_reads(read1, f_out.name, unmapped_reads)
        _filter_reads(read2, f_out.name, unmapped_reads)
    
    readprocessing2 = [
        f"bwa mem {REF} {dual_read1} {dual_read2} > {readdir}/target_{tm.edgeid}.sam",
        f"{samtools} view -bT {REF} {readdir}/target_{tm.edgeid}.sam > {readdir}/target_{tm.edgeid}.bam",
        f"{samtools} sort {readdir}/target_{tm.edgeid}.bam > {readdir}/target_{tm.edgeid}_sorted.bam",
        f'{samtools} view -F 4 -f 8 {readdir}/target_{tm.edgeid}.bam |' + " awk '{print $1, $2, $4, $8, $5}' > lmapped.txt",
        f'{samtools} view -F 8 -f 4 {readdir}/target_{tm.edgeid}.bam |'  " awk '{print $1, $2, $4, $8, $5}' > rmapped.txt",
        f'{samtools} view -f 0x2 {readdir}/target_{tm.edgeid}.bam |' + " awk '{print $1, $2, $4, $8, $5}' > fullmapped.txt"]


    for line in readprocessing2:
        subprocess.call(line, shell = True, cwd=readdir)

    data = f"{readdir}/fullmapped.txt"
    df = pd.read_csv(data, sep = " ", header = None, 
                    names=["read_name", "flag", "pos", "mpos", "mapq"])

    df["flag"] = df["flag"].apply(lambda x: bin(x)[2:].zfill(12))
    cols = df["flag"].apply(lambda x: list(x[-11:]))

    bnames = ["paired", "proper", "unmapped_q", "mate_unmapped", "strand", "mate_strand",
            "first_read", "second_read", "not_prim", "QC", "op_PCR"][::-1]
    df = pd.concat([df, pd.DataFrame(cols.tolist(), columns = bnames).astype(int)], axis = 1)
    df.drop(['flag', 'QC', 'op_PCR', 'not_prim'],
            axis = 1, inplace = True)

    df = df[df["proper"] == 1] 
    df = df[df["first_read"] == 1]
    # df.dropna(subset=["l", "r"], inplace=True) #FIXME dangerous
    df[["l", "r"]] = df.apply(compute_lr, axis = 1)

    threshold = 500 # check
    sus_df = df[(df["l"] > df["r"]) | (df["r"] - df["l"] > threshold)]

    if log != "":
        with open(log, "a+") as f:
            f.write(f"Generated {d['numreads']} reads")

    return df, sus_df


def readstep2(d, readdir, rp, tm, aug_fp, verbose = False, log = ""):
    cprint = lambda x: print(x) if verbose else ""
    cellname = tm.edgeid.split(".")[0]

    read1 = f"{readdir}/{rp}.bwa.read1.fastq"
    read2 = f"{readdir}/{rp}.bwa.read2.fastq"

    line =f"dwgsim -N {d['numreads']} -1 {d['readlen']} -2 {d['readlen']} -d {d['outerdist']} -S {d['readmode']} {aug_fp} {rp}"
    subprocess.call(line, shell = True, cwd=readdir)

    compute_lr = _compute_lr(d['readmode'])

    # Define input and output files
    cwd = readdir
    read1 = f"{readdir}/{rp}.bwa.read1.fastq.gz"
    read2 = f"{readdir}/{rp}.bwa.read2.fastq.gz"
    dual_read1 = f"{readdir}/filtered_read1.fastq"
    dual_read2 = f"{readdir}/filtered_read2.fastq"
    unmapped = f"{readdir}/unmapped.fastq"

    # 1. Extract read names from both files
    cprint("Extracting read names...")
    read1_names = _get_read_names(read1, cellname)
    read2_names = _get_read_names(read2, cellname)

    # 2. Identify dually mapped and unmapped reads
    dually_mapped = read1_names & read2_names
    unmapped_reads = (read1_names | read2_names) - dually_mapped

    cprint(f"Found {len(dually_mapped)} dually mapped reads.")
    cprint(f"Found {len(unmapped_reads)} unmapped reads.")

    # 4. Filter dually mapped reads
    cprint("Filtering dually mapped reads...")
    _filter_reads(read1, dual_read1, dually_mapped)
    _filter_reads(read2, dual_read2, dually_mapped)

    # 5. Filter unmapped reads (combine both ends)
    cprint("Filtering unmapped reads...")
    with open(unmapped, 'w') as f_out:
        _filter_reads(read1, f_out.name, unmapped_reads)
        _filter_reads(read2, f_out.name, unmapped_reads)
    
    readprocessing2 = [
        f"bwa mem {REF} {dual_read1} {dual_read2} > {readdir}/target_{tm.edgeid}.sam",
        f"{samtools} view -bT {REF} {readdir}/target_{tm.edgeid}.sam > {readdir}/target_{tm.edgeid}.bam",
        f"{samtools} sort {readdir}/target_{tm.edgeid}.bam > {readdir}/target_{tm.edgeid}_sorted.bam",
        f'{samtools} view -F 4 -f 8 {readdir}/target_{tm.edgeid}.bam |' + " awk '{print $1, $2, $4, $8, $5}' > lmapped.txt",
        f'{samtools} view -F 8 -f 4 {readdir}/target_{tm.edgeid}.bam |'  " awk '{print $1, $2, $4, $8, $5}' > rmapped.txt",
        f'{samtools} view -f 0x2 {readdir}/target_{tm.edgeid}.bam |' + " awk '{print $1, $2, $4, $8, $5}' > fullmapped.txt"]


    for line in readprocessing2:
        subprocess.call(line, shell = True, cwd=readdir)

    
    threshold = 500
    compute_lr = _compute_lr(1)
    data = f"{readdir}/fullmapped.txt"

    rows = []
    with open(data, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            read_name = parts[0]
            flag = int(parts[1])
            pos = int(parts[2])
            mpos = int(parts[3])
            mapq = int(parts[4])
            rows.append((read_name, flag, pos, mpos, mapq))

    bnames = ["paired", "proper", "unmapped_q", "mate_unmapped", "strand", "mate_strand",
            "first_read", "second_read", "not_prim", "QC", "op_PCR"][::-1]
    flags_bin  = [bin(flag)[2:].zfill(12)[-11:] for _, flag, _, _, _ in rows]

    decoded = [{name: int(bit) for name, bit in zip(bnames, bits)}  for bits in flags_bin]

    filtered_rows = [
        (*row, decoded[i])
        for i, row in enumerate(rows)
        if decoded[i]["proper"] == 1 and decoded[i]["first_read"] == 1]
    
    # ----- Compute l, r -----
    results = []
    for row in filtered_rows:
        read_name, flag, pos, mpos, mapq, flag_dict = row
        l_r = compute_lr({
            "pos": pos,
            "mpos": mpos,
            "strand": flag_dict["strand"]
        })
        l = l_r["l"]
        r = l_r["r"]
        results.append((read_name, l, r, mapq)) 

    # ----- Apply suspect filter -----
    suspects = [
        (read_name, l, r, mapq)
        for read_name, l, r, mapq in results
        if (l > r or (r - l > threshold))
    ]

    if log != "":
        with open(log, "a+") as f:
            f.write(f"Generated {len(filtered_rows)} mapped reads, found {len(suspects)} jump events.")

    return results, suspects
    

