import numpy as np
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots 
import plotly.graph_objs as go
import matplotlib.pyplot as plt

def visualize_coverage(coverage, length):
        """length: length of reference (index)"""
        counts = coverage
        total_count = sum(counts)
        fig = make_subplots(rows = len(counts) + 1, cols = 1, 
                            shared_xaxes = True, vertical_spacing = 0.1)
        for i in range(len(counts)):
            fig.add_trace(go.Scatter(x = np.arange(length), y = counts[i]), 
                          row = i + 1, col = 1)
            fig.update_xaxes(title_text = f'Population {i + 1} Coverage', row = i + 1, col = 1)
        fig.add_trace(go.Scatter(x = np.arange(length), y = total_count),
                      row = len(counts) + 1, col = 1)
        fig.update_xaxes(title_text = "Total Coverage", row = len(counts) + 1, col = 1)
        fig.update_layout(height=600, width=600, title_text="Coverage Representation")
        fig.show()

def visualize_coords(coords):
    traces = []

    # downsample data
    sample_range = min([len(coord) for coord in coords])
    sample_index = np.random.choice(sample_range, size = 1000, replace = False)
    coords = [[coord[i] for i in sample_index] for coord in coords]

    # plot scatter plot
    for i in range(len(coords)):
        traces.append(go.Scatter(x=[x[0] for x in coords[i]], y=[x[1] for x in coords[i]], 
                                    name = f'Population {i}', mode = 'markers', marker = dict(size = 4, opacity = 0.5)))
        
    layout = go.Layout(title='Paired-End Representation', xaxis=dict(title='left read starting position'), yaxis=dict(title='right read starting position'))
    fig = go.Figure(data = traces, layout = layout)
    fig.show()

def save_total_counts(filename, total):
    pd.DataFrame({"counts": total}).to_csv(filename)

def save_coords(filename, coords):
     pd.DataFrame(sum(coords, []), columns = ['left read', 'right read']).to_csv(filename)