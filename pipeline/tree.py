import numpy as np 
from dataclasses import dataclass, field
from typing import Optional, List
import copy 

@dataclass
class Cell:
    # node properties
    cellname: str 
    children: List["Cell"] = field(default_factory=list)
    mutations: List["Mutation"] = field(default_factory=list)
    parent: Optional["Cell"] = None
    founder: Optional["Cell"] = None

    # sequence interests

    # not implemented
    p: float = None #FIXME relative proportion in mixture, when relevant

    def parent_path(self):
        parent = self.parent 
        path = []
        while parent is not None:
            path.append(parent)
            parent = parent.parent 
        
        return path[::-1]
    
    def all_mutations(self, allele):
        path = self.parent_path()
        ms = []
        for p in path:
            for m in p.mutations:
                if m.allele == allele:
                    ms.append(m)

        for m in self.mutations:
            if m.allele == allele:
                ms.append(m)

        return ms

    def sequence(self, allele, out=""):
        if self.founder is not None:
            seq = copy.deepcopy(self.founder.ref)
        else:
            seq = copy.deepcopy(self.ref)
        parent_path = self.parent_path()
        if len(parent_path) > 0:
            for p in parent_path:
                for m in p.mutations:
                    if m.allele == allele:
                        seq = m.mutate_relative(seq)
        # self mutations
        for m in self.mutations:
            if m.allele == allele:
                seq = m.mutate_relative(seq)
        
        return seq 

    def __str__(self):
        return self.cellname

@dataclass
class Mutation:
    edgeid: int 
    allele: int 
    num_copies: int
    relstart: int 
    refstart: int 
    length: int 
    d: dict # passing in cnasim utils, etc. 

    def __post_init__(self):
        self.left = self.relstart
        self.right = self.relstart + self.length
        
        self.l = self.relstart // self.d["resolution"]
        self.r = self.l + self.length // self.d["resolution"]
    
    def mutate_relative(self, seq):
        return seq[:self.l] + self.num_copies * seq[self.l:self.r] + seq[self.r:]



@dataclass
class Mutation2:
    #FIXME only one mutation supported per generation
    edgeid: int
    start: Cell 
    end: Cell 
    allele: int
    num_copies: int 
    relstart: int
    refstart: int
    length: int
    d: dict # passing cnasim arguments, etc. 
    # FIXME some absolute information

    def __post_init__(self):
        self.left = self.relstart
        self.right = self.relstart + self.length

        # provide resolution 
        self.l = self.relstart // self.d["block_length"]
        self.r = self.l + self.length // self.d["block_length"]

        def _mutate_relative(self):
            def mutate(seq):
                if self.num_copies == 0:
                    return seq[:self.l] + seq[self.r]
                else:
                    region = seq[self.l:self.r]
                    return seq[:self.l] + self.num_copies * region + seq[self.l]
            return mutate 
        
        self.mutate = _mutate_relative
        
        # assumption - tree traversal has a single route to this edge
        def _mutate_absolute(self):
            transforms = [self.mutate]
            cnode = self.start 
            while cnode.parent is not None:
                transforms.append(cnode.parent_path.mutate)
                cnode = cnode.parent      
            for t in transforms[::-1]: 
                mut = lambda seq: t(seq)
            return mut
        self.absolute_mutation = _mutate_absolute
        
    def __str__():
        pass

def parse_newick_structure(nwk):
    #NB discarding distance info on build
    stack = []
    root = []  # This will store the final parsed structure
    current = root  # This is the current level we're working on

    token = ""
    for char in nwk:
        if char == "(":
            # Start a new subtree
            subtree = []
            current.append(subtree)
            stack.append(current)  # Save the current level
            current = subtree  # Move to the new level
        elif char == ",":
            if token:
                nodename = token.strip().split(":")[0]
                current.append(nodename)  # Store the node name
                token = ""  # Reset for the next node
        elif char == ")":
            if token:
                nodename = token.strip().split(":")[0]
                current.append(nodename)  # Store the last node before closing the subtree
                token = ""  # Reset for the next node
            current = stack.pop()  # Pop the parent level from the stack
        elif char == ";":
            if token:
                nodename = token.strip().split(":")[0]
                current.append(nodename)  # Store the last node before end
        else:
            token += char  # Collect node name characters
    


    return root

# Test the function
# nwk = "(cell1:0.7272168501672318,(cell2:0.0725908101359075,(cell3:0.04753208917522399,(cell4:0.04267063059524195,cell5:0.04267063059524195)ancestor3:0.00486145857998203)ancestor2:0.02505872096068351)ancestor1:0.6546260400313243)founder;"
# parsed_structure = parse_newick_structure(nwk)

def parse_focal_subset(focal):
    #FIXME accessible class
    return focal.split()



@dataclass
class Population:
    #FIXME flexibility with initializing your own 
    focal: dict
    ref: str
    nwk_str: str
    root: Cell

    def __post_init__(self):

        self.parsed = parse_newick_structure(self.nwk_str)

        # FIXME shallow copying
        parsed = self.parsed[::]

        traversal = []
        parent = None

        # pass one - unencoded-mutation case
        while len(self.parsed) > 0:
            cname = parsed.pop().split[':'][0]


            if cname == "founder":
                founder = Cell(cname)
                founder.ref = self.ref
                traversal.append[founder]
                
            else:
                cnode = Cell(cname, founder = founder, parent = parent)
                traversal.append(cnode)
                parent.children.append(cnode)

            if type(parsed[-1]) == list:
                parent = cnode

        # pass two - focal encoding plus edge connections
        for i in range(len(traversal)):
            #BUG assumption of one event unreasonable
            ordered_events = parse_focal_subset()
    
    def __str__(self):
        print(":/") 

def build_tree(traversal, focal, d):
    """
    Build a tree from a right-to-left traversal list.
    
    At each level the list is assumed to consist of one or more groups encoded as:
    
       [ immediate_children, parent ]
    
    where immediate_children is exactly the element immediately to the left of the parent's name.
    (If that element is a list of strings, it yields leaf nodes; if it is a mixed list,
    it is processed recursively in the same way.) Any elements further left form another group
    at the same level.
    
    The root is given as the last element of the overall traversal.
    """
    # Create the root (founder)
    root_name = traversal[-1]
    root = Cell(cellname=root_name)

    celldict = {root_name: root}
    
    def _populate_mutations(node, focal):
        _id = 0
        for row in focal.get(node.cellname, []):
            allele = int(row[0])
            tumor_start = int(row[1]) - 1
            ref_start = int(row[2]) - 1
            length = int(row[3])
            copies = int(row[4] == "gain") * (1 + int(row[5]))
            node.mutations.append(Mutation(
                edgeid=f'{node.cellname}.{_id}',
                allele=allele,
                num_copies=copies,
                relstart=tumor_start,
                refstart=ref_start,
                length=length,
                d=d))
            _id += 1
    _populate_mutations(root, focal)
    def _process_group(lst, parent):
        """
        Process a list 'lst' assumed to encode one or more groups.
        Returns nothing (nodes are attached to parent).
        """
        i = len(lst) - 1
        while i >= 0:
            # The parent's name is always at the end of a group.
            # If we have at least two elements, then lst[i] is the parent's name
            # and lst[i-1] is its immediate children.
            if i - 1 >= 0:
                group_parent = lst[i]  # Expected to be a string
                children_elem = lst[i-1]
                # Create the parent node for this group:
                parent_node = Cell(cellname=group_parent, parent=parent, founder=root)
                _populate_mutations(parent_node, focal)
                parent.children.append(parent_node)
                celldict[group_parent] = parent_node
                # Now, process the immediate children.
                if isinstance(children_elem, list):
                    # If it's a pure list of strings, create leaves.
                    if all(isinstance(x, str) for x in children_elem):
                        for child_name in children_elem:
                            child = Cell(cellname=child_name, parent=parent_node, founder=root)
                            _populate_mutations(child, focal)
                            parent_node.children.append(child)
                            celldict[child_name] = child
                    else:
                        # Mixed list: process recursively.
                        _process_group(children_elem, parent_node)
                elif isinstance(children_elem, str):
                    # A single leaf.
                    child = Cell(cellname=children_elem, parent=parent_node, founder=root)
                    _populate_mutations(child, focal)
                    parent_node.children.append(child)
                    celldict[children_elem] = child
                else:
                    raise ValueError("Unexpected type in group (expected list or string).")
                i -= 2  # We've consumed this group.
            else:
                # If there's a leftover element that doesn't form a complete group, treat it as a leaf.
                elem = lst[i]
                if isinstance(elem, str):
                    leaf = Cell(cellname=elem, parent=parent, founder=root)
                    _populate_mutations(leaf, focal)
                    parent.children.append(leaf)
                    celldict[elem] = leaf
                elif isinstance(elem, list):
                    # Process the list as a pure group of leaves.
                    for s in elem:
                        leaf = Cell(cellname=s, parent=parent, founder=root)
                        _populate_mutations(leaf, focal)
                        parent.children.append(leaf)
                        celldict[elem] = leaf
                i -= 1

    def _recursive_build(subtree, parent):
        """
        Process a subtree level.
        If the subtree is a pure list of strings, attach them as leaves.
        Otherwise, treat the subtree as one or more groups.
        """
        if isinstance(subtree, list) and all(isinstance(x, str) for x in subtree):
            # Pure list of leaves.
            for s in subtree:
                leaf = Cell(cellname=s, parent=parent, founder=root)
                _populate_mutations(leaf, focal)
                parent.children.append(leaf)
        elif isinstance(subtree, list):
            _process_group(subtree, parent)
        else:
            # Not a list: it should be a string.
            leaf = Cell(cellname=subtree, parent=parent, founder=root)
            _populate_mutations(leaf, focal)
            parent.children.append(leaf)
    
    # The main traversal's first element is the top-level subtree.
    _recursive_build(traversal[0], root)
    return root, celldict

# Function to print tree structure
def print_tree(node, level=0):
    print("  " * level + node.cellname)
    for child in node.children:
        print_tree(child, level + 1)

def check_non_overlapping(cell, allele):
    ms = cell.all_mutations(allele)

    if len(ms) <= 1:
        return ms
    
    # Create intervals (start, end) and keys (i, i)
    vals = [(m.relstart, m.relstart + m.num_copies * m.length) for m in ms]
    parents = [m.edgeid.split('.')[0] for m in ms]
    keys = [(i, i) for i in range(len(ms))]
    
    # Zip and flatten by sorting on interval points
    zipped = list(zip(vals, keys))
    flat = sorted([(point, i) for (start, end), (i, _) in zipped for point in (start, end)])
    
    # Identify non-overlapping indices
    non_overlapping = set()
    last_id = None
    for _, i in flat:
        if i == last_id:
            non_overlapping.add(i)  # If i appears back-to-back, it's non-overlapping
        last_id = i
    
    return [ms[i] for i in non_overlapping]  # Return sorted list of non-overlapping indices
