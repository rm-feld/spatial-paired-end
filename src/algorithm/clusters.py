import pandas as pd
from ..seqs_and_reads.consts import *
import sklearn 
import matplotlib.pyplot as plt
import numpy as np

# filter out 'main region' 
def only_mutations(df):
    return df[((df["left read"] + 2 * READ_LEN + READ_SEP[0]) > df["right read"])|((df["left read"] + 2 * READ_LEN + READ_SEP[1]) < df["right read"])]


def square_density(df, centroid, sidelength):
    step = sidelength // 2
    x, y = centroid
    rdf = df[(x - step < df["left read"]) & (df["left read"] < x + step) & (y - step < df["right read"]) & (df["right read"] < y + step)]
    return rdf.shape[0] / (step * 2) ** 2

def circle_density(df, centroid, radius):
    pass

def prop_calculation(dena, denb, denc, mtype, baseline):
    if mtype == "deletion":
        pass
    else:
        pass

