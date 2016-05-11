import pandas as pd
import numpy as np
import json
import os

data_root = '/home/cbarnes/work/code/connectome/paper/metrics/graphs/di_layers'
out_root = '/home/cbarnes/work/code/connectome/paper/metrics/csvs'

metrics = [
    'assortativity',
    'density',
    'global_efficiency',
    'mean_betweenness_centrality',
    'mean_clustering',
    'modularity',
    'path_length',
    'small_worldness',
    'transitivity',
    'weighted_mean_clustering'
]

graph_names = [
    ('syn', 'Syn'),
    ('gj', 'Gap'),
    ('gj-syn', 'Wired'),
    ('ma', 'MA'),
    ('gj-ma-syn', 'All'),
    ('np', 'NP'),
    ('gj-ma-np-syn', 'All w/ NP'),
    ('gj-ma-syn2', 'All (real wired + null MA)')
]

cols = ['real'] + ['null_{:02.0f}'.format(i) for i in range(100)]


def filename_iter(nreps=np.inf, ext='.npy'):
    i = 0
    while i < nreps:
        yield '{:03}{}'.format(i, ext)
        i += 1


def from_json(path, key):
    with open(path) as f:
        return json.load(f)[key]


for metric in metrics:
    rows = []
    for graph_name, index_name in graph_names:
        row = [from_json(os.path.join(data_root, graph_name, 'adj.json'), metric)]

        for filename in filename_iter(100, '.json'):
            path = os.path.join(data_root, graph_name, 'controls', filename)
            row.append(from_json(path, metric))

        rows.append(row)

    df = pd.DataFrame(data=rows, index=[tup[1] for tup in graph_names], columns=cols)
    df.to_csv(os.path.join(out_root, metric + '.csv'))

with open(os.path.join(out_root, 'meta.txt'), 'w') as f:
    f.write('''\
assortativity: full-full degree
density:
global_efficiency:
mean_betweenness_centrality: unweighted
mean_clustering: mean nodewise clustering coefficient, NOT global transitivity
modularity: maximum modularity - probably don't trust this as bctpy gives different answers to BCT due to differences in underlying linear algebra implementation
path_length:
small_worldness: probably don't trust this as
transitivity: use this as clustering
weighted_mean_clustering: mean nodewise clustering coefficient weighted by full degree of nodes
''')

print('done')