import numpy as np
import bct
import os
from scipy import stats
try:
    from metrics.file_tools import filename_iter
    from metrics.shared import set_seeds, push_exceptions
except (ImportError, SystemError):
    from file_tools import filename_iter
    from shared import set_seeds, push_exceptions

set_seeds(1)

DATA_ROOT = 'graphs/di_layers/'

real_phys = np.load(os.path.join(DATA_ROOT, 'gj-syn', 'adj.npy'))
real_ma = np.load(os.path.join(DATA_ROOT, 'ma', 'adj.npy'))

struc, Q = bct.modularity_dir(real_phys)


def adj_to_edgelist(adj):
    return np.argwhere(adj)


def edge_in_module(edge, struc):
    src, tgt = edge
    return struc[src] == struc[tgt]


def edges_in_module(edges, struc, proportion=True):
    in_mod = [edge_in_module(edge, struc) for edge in edges]
    if proportion:
        return sum(in_mod) / len(edges)
    else:
        return sum(in_mod)


def random_edges(n):
    edge_set = set()
    while len(edge_set) < n:
        src, tgt = np.random.randint(0, 302, 2)
        if src == tgt:
            continue
        edge_set.add((src, tgt))

    return sorted(edge_set)


def get_real_prop():
    return edges_in_module(adj_to_edgelist(real_ma), struc)


def get_random_prop(std=False, reps=100):
    vals = [edges_in_module(random_edges(real_ma.sum()), struc) for _ in range(reps)]

    return (np.mean(vals), np.std(vals)) if std else np.mean(vals)


def null_iter():
    for filename in filename_iter(100):
        yield np.load(os.path.join(DATA_ROOT, 'ma', 'controls', filename))


def get_null_prop(std=False):
    vals = [edges_in_module(adj_to_edgelist(null_graph), struc) for null_graph in null_iter()]
    return (np.mean(vals), np.std(vals)) if std else np.mean(vals)


def get_zscores(data, absolute=False):
    reals_controls = zip([d['real'] for _, d in data], [d['control'] for _, d in data])
    if absolute:
        return [abs((real - np.mean(controls))/np.std(controls)) for real, controls in reals_controls]
    else:
        return [(real - np.mean(controls))/np.std(controls) for real, controls in reals_controls]


def data_to_pval_strs(data):
    pvals = stats.norm.sf([abs(val) for val in get_zscores(data)])*2  # 2-sided
    return ['p={:.2e}'.format(val) for val in pvals]


if __name__ == "__main__":
    real_prop = get_real_prop()
    rand_prop, rand_std = get_random_prop(std=True)
    rand_pval = (1 - stats.norm.cdf(real_prop, loc=rand_prop, scale=rand_std))
    null_prop, null_std = get_null_prop(std=True)
    null_pval = (1 - stats.norm.cdf(real_prop, loc=null_prop, scale=null_std))

    print('''
real_prop: {real_prop}
rand_prop: {rand_prop}, p = {rand_pval}
null_prop: {null_prop}, p = {null_pval}
'''.format(**locals()))

    print('There are {} modules in the phys network'.format(len(np.unique(struc))))
