import bct
import numpy as np
import os
import json
import networkx as nx
import warnings
import random
try:
    from metrics.file_tools import filename_iter
    from metrics.shared import set_seeds, push_exceptions
except (ImportError, SystemError):
    from file_tools import filename_iter
    from shared import set_seeds, push_exceptions

set_seeds()

N_CONTROLS = 100
EXTRASYN_TP = 'ma'
DATA_ROOT = 'graphs'


def gen_controls(real_phys, extrasyn_root, n_controls=N_CONTROLS):
    """
    Get controls from files (randomised by edge swap)
    """
    for filename in filename_iter(n_controls):
        path = os.path.join(extrasyn_root, 'controls', filename)
        try:
            control_extrasyn = np.load(path)
        except FileNotFoundError:
            warnings.warn(
                'Tried to get modularity for {} controls but only found up to {}'.format(n_controls, filename)
            )
            raise StopIteration
        yield bct.binarize(real_phys + control_extrasyn)


def actual_random(edges, nodes=302):
    """
    Generate completely random undirected non-self-looping graph with given number of edges
    """
    W = np.zeros((nodes, nodes))
    while W.sum() < edges:
        pair = sorted(set(np.random.randint(0, 302, 2)))
        if len(pair) > 1 and W[pair[0], pair[1]] == 0:
            W[pair[0], pair[1]] = 1
    return W + W.T


def gen_rands(real_phys, extrasyn_root, n_controls=N_CONTROLS):
    cont = np.load(os.path.join(extrasyn_root, 'adj.npy'))
    edges = np.sum(np.triu(cont, 1))
    for _ in range(n_controls):
        yield bct.binarize(real_phys + actual_random(edges))


def shuffle_nodes(adj):
    """
    Generate random graph by shuffling the node labels
    """
    arr = adj.copy()
    g = nx.from_numpy_matrix(arr)
    nodelist = g.nodes()
    random.shuffle(nodelist)
    return np.array(nx.to_numpy_matrix(g, nodelist))


def gen_shuffled(real_phys, extrasyn_root, n_controls=N_CONTROLS):
    cont = np.load(os.path.join(extrasyn_root, 'adj.npy'))
    for _ in range(n_controls):
        yield bct.binarize(real_phys + shuffle_nodes(cont))


@push_exceptions
def main():
    print('MODULARITY CHANGE')
    for source in ['ac', 'ww']:
        print('  ' + source)
        for ma_tp in ['str', 'wk']:
            print('    ' + ma_tp)
            extrasyn_tp = EXTRASYN_TP
            if extrasyn_tp == 'ma':
                extrasyn_tp += '_' + ma_tp

            extrasyn_root = os.path.join(DATA_ROOT, extrasyn_tp)

            real_phys = np.load(os.path.join(DATA_ROOT, 'gj-syn_{}'.format(source), 'adj.npy'))
            struc, mod = bct.modularity_und(real_phys)

            # generate real values
            real_extrasyn = np.load(os.path.join(extrasyn_root, 'adj.npy'))
            real_comb = bct.binarize(real_phys + real_extrasyn)
            real_comb_struc, real_comb_mod = bct.modularity_und(real_comb, kci=struc)
            assert np.allclose(real_comb_struc, struc)

            d = {
                'real_phys': mod,
                'real_comb': real_comb_mod,
            }
            # loop through the randomisation algorithms, generating the modularities for each
            name_alg = {'control': gen_controls,
                        'random': gen_rands,
                        'shuffled': gen_shuffled}
            for name, alg in sorted(name_alg.items()):
                print('      Randomising by ' + name)
                mods = []
                for control in alg(real_phys, extrasyn_root):
                    print('.', end='')
                    comb_struc, comb_mod = bct.modularity_und(control, kci=struc)
                    assert np.allclose(comb_struc, struc)
                    mods.append(comb_mod)
                d2 = {'{}_comb'.format(name): mods,
                      '{}_comb_mean'.format(name): np.mean(mods),
                      '{}_comb_std'.format(name): np.std(mods)}
                d.update(d2)
                if name == 'shuffled':
                    1+1
                print()

            # write to file
            out_root = 'mod_plots/{}_{}'.format(source, extrasyn_tp)
            os.makedirs(out_root, exist_ok=True)
            with open(os.path.join(out_root, 'data.json'), 'w') as f:
                json.dump(d, f, indent=2, sort_keys=True)

if __name__ == '__main__':
    main()
