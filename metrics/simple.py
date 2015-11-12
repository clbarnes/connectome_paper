"""
WW only. No neuropeptides. No controls. No weak MA edges. Physical vs Combined. Undirected, unweighted only.

density TICK
path length TICK
global efficiency TICK
clustering TICK
modularity TICK
small-worldness REQUIRES RANDOMS
degree distribution TICK
assortivity TICK
centrality TICK
"""

from multiplex import MultiplexConnectome
import connectome_utils as utl
import numpy as np
import networkx as nx
import bct
from collections import Counter
import json
import scipy.io

whole_path = '/home/cbarnes/work/code/connectome/paper/real_graphs/whole_no_np/strong_only/ww/graph.json'
whole_graph = utl.json_deserialise(whole_path)
whole_C = MultiplexConnectome(whole_graph)
nodelist = sorted(whole_C.whole.nodes())


def networkx_to_und_unwe_mat(G):
    return np.array(nx.to_numpy_matrix(
        G.to_undirected(),
        nodelist,
        weight=None,
        multigraph_weight=min
    ))


def networkx_to_dir_unwe_mat(G):
    return np.array(nx.to_numpy_matrix(
        G,
        nodelist,
        weight=None,
        multigraph_weight=min
    ))


def make_phys_comb(multcon=whole_C):
    return networkx_to_und_unwe_mat(multcon.compose('GapJunction', 'Synapse')), networkx_to_und_unwe_mat(multcon.whole)


def density(adj):
    return bct.density_und(adj)[0]


def path_length(adj):
    return bct.charpath(bct.distance_bin(adj))[0]


def global_efficiency(adj):
    return bct.charpath(bct.distance_bin(adj))[1]


def clustering(adj):
    return bct.clustering_coef_bu(adj).mean()


def transitivity(adj):
    return bct.transitivity_bu(adj)


def modularity(adj):
    return bct.modularity_und(adj)[1]


def assortativity(adj):
    return bct.assortativity_bin(adj, 0)


def betweenness_centrality(adj):
    return bct.betweenness_bin(adj).mean()


def degree_dist(adj):
    """
    There are y nodes with a degree of x
    """

    counts = Counter(np.array(adj).sum(axis=1))
    xy_pairs = []
    for x, y in sorted(counts.items()):
        xy_pairs.append([x, y])

    return xy_pairs

algorithms = {
    'density': density,
    'path length': path_length,
    'global efficiency': global_efficiency,
    'clustering coefficient': clustering,
    'transitivity': transitivity,
    'maximum_modularity': modularity,
    'assortativity': assortativity,
    'betweenness centrality': betweenness_centrality,
    'degree distribution': degree_dist
}

# phys, comb = make_phys_comb()

results = dict()

phys = scipy.io.loadmat('/home/cbarnes/work/code/bctpy_test/data/phys_ww_mats.mat')['und']

np.savetxt('example.txt', phys, delimiter=',', fmt='%i')

for name, algorithm in algorithms.items():
    results[name] = dict()
    results[name]['physical'] = algorithm(phys)
    # results[name]['combined'] = algorithm(comb)

with open('simple_results.json', 'w') as f:
    json.dump(results, f, sort_keys=True, indent=2)

print('ready')