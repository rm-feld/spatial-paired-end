import numpy as np
import random

class Node:
    def __init__(self, name, val = None, parent = None):
        self.name = name 
        self.val = val
        self.parent = parent
        self.children = []

class Tree:
    def __init__(self, nodelist, parent):
        self.parent = parent 
        self.nodelist = nodelist
    
    def __str__(self):
        #FIXME implement vals = True
        output = "---------\n"
        for node in self.nodelist:
            val = ", vals: " + str(node.val) if node.val != None else ""
            output += f"name: {node.name}{val}\n"
            parent = "root" if node.parent == None else node.parent.name
            output += f"\tparent: {parent}\n"
            output += f"\tchildren: {node.children}\n"
        output += "---------"
        return output
    
    def get_node(self, i):
        return self.nodelist[i]
    
    def get_node_ancestry(self, i):
        output = []
        cnode = self.nodelist[i]
        while cnode.parent is not None:
            output.append(cnode.parent.name)
            cnode = cnode.parent
        return output
    
    def matrix_form(self):
        n = len(self.nodelist)
        A = np.eye(n)
        for i in range(n):
            ancestors = self.get_node_ancestry(i)
            for j in ancestors:
                if j != 0:
                    A[i][j] = 1
        A[0][0] = 0
        return A

def build_random_tree(n, vals = None):
    nodes = [Node(0)]
    for i in range(1, n):
        p_id = random.randint(0, len(nodes) - 1)
        node_i = Node(i, parent = nodes[p_id])
        nodes[p_id].children.append(i)
        nodes.append(node_i)
    
    return Tree(nodes, parent = nodes[0])

class Graph:
    def __init__(self, nodearr, edgearr):
        self.nodearr = nodearr
        self.edgearr = edgearr

class Edge:
    def __init__(self, n1, n2):
        self.n1 = n1 
        self.n2 = n2
        self.avoid = False

class GNode:
    def __init__(self, p, i):
        self.p = p
        self.i = i 

class Seq:
    def __init__(self, rand = True, n = None, m = None, p = None, tree = None):
        if rand == True:
            self.tree = build_random_tree(n)
            self.m = np.concatenate((np.asarray([1]), np.random.randint(2, 5, size = n - 1)), axis = None)
            self.p = np.random.dirichlet(np.ones(n)) + 1/(n ** 2) # 1/n^2 is a normalizing term (#TODO: write out/determine current viability of small case)
            self.p = self.p / np.sum(self.p)
            self.n = n
        else:
            self.tree = tree
            self.m = m 
            self.p = p 
            self.n = len(m)


        self.I = np.multiply(self.tree.matrix_form(), self.m)
        self.I[self.I == 0] = 1
        self.X = (self.I.T * self.p).T
        self.D = np.asarray([np.sum(self.X[:, i]) - 1 for i in range(n)])


    
    def convert_to_graph(self):
        nodearr = []
        edgearr = []

    
     
# seq = Seq(n = 6)

# print(seq.tree)

# print(seq.I)
# print("p:")
# print(seq.p)
# print(seq.X)
# print()
# print()

# print(seq.D)
