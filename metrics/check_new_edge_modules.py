import numpy as np
import bct
import os
try:
    from metrics.file_tools import filename_iter
    from metrics.shared import set_seeds, push_exceptions
except (ImportError, SystemError):
    from file_tools import filename_iter
    from shared import set_seeds, push_exceptions

DATA_ROOT = 'graphs/di_layers/'

real_phys = np.load(os.path.join(DATA_ROOT, 'gj-syn', 'adj.npy'))

struc, Q = bct.modularity_dir(real_phys)


def edge_in_module(edge, struc):
    src, tgt = edge
    return struc[src] == struc[tgt]


def edges_in_module(edges, struc, proportion=True):
    in_mod = [edge_in_module(edge, struc) for edge in edges]
    if proportion:
        return sum(in_mod) / len(edges)
    else:
        return sum(in_mod)


# todo: compute with completely random edges, and with null MA graphs, and with real MA graph