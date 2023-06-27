import numpy as np
from .consts import *
from .mutate import random_mutations
from .sequencing import *

def generate_base_seq():
    l = np.random.randint(FRAGMENT_LENGTH[0], FRAGMENT_LENGTH[1])
    return np.arange(l)

def generate_pops():
    pass

class Populations(object):
    def __init__(self, is_random = True):
        if not is_random: 
            self.pops = 0
            self.pop_mix = []
        else:
            self.pops = 4
            self.pop_mix = [0.3, 0.2, 0.2, 0.3]
    def generate(self):
        self.ref = generate_base_seq()
        pop_data = [self.ref]
        mutation_data = []
        read_data = []
        coord_data = []
        coverage_data = []

        # generate sequences
        for i in range(self.pops - 1):
            new_string, mutations = random_mutations(pop_data[i])
            pop_data.append(new_string)
            mutation_data.append((i + 1, mutations))
        
        # generate reads and coverage
        for i in range(self.pops):
            i_coverage = int(self.pop_mix[i] * COVERAGE)
            reads, coords = paired_end_reads(pop_data[i], i_coverage)
            coverage = find_coverage(reads, self.ref)
            read_data.append(reads)
            coord_data.append(coords)
            coverage_data.append(coverage)

        # organize output
        dic = {
            "strands": pop_data,
            "mutations": mutation_data,
            "reads": read_data, 
            "coordinates": coord_data, 
            "coverage": coverage_data,
            "total": np.sum(coverage_data, axis = 0)
        }

        self.dic = dic 
        return dic
    
    def add_seq(self, ref, seq, mutations, pop_mix):
        #TODO: let us add more than one sequence
        if self.pops == 0:
            self.pops = 2 
            self.pop_mix = pop_mix
            read_data = [],
            coord_data = [], 
            coverage_data = []
            pop_data = [ref, seq]
            mutation_data = [(1, mutation) for mutation in mutations]
            for i in range(self.pops):
                i_coverage = int(self.pop_mix[i] * COVERAGE)
                reads, coords = paired_end_reads(self.pop_data[i], i_coverage)
                coverage = find_coverage(reads)
                read_data.append(reads)
                coord_data.append(coords)
                coverage_data.append(coverage)
            # organize output
            dic = {
                "strands": pop_data,
                "mutations": mutation_data,
                "reads": read_data,
                "coordinates": coord_data,
                "coverage": coverage_data,
                "total": np.sum(coverage_data, axis = 0)}
            self.dic = dic 
 

        



    


