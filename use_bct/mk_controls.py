import connectome_utils as utl
from multiplex import MultiplexConnectome
import bct
import os
import numpy as np
import networkx as nx
import shutil

reps = 100
swap_prop = 10

ac_comb = utl.json_deserialise('/home/cbarnes/work/code/connectome/paper/real_graphs/whole_no_np/strong_only/ac/graph'
                               '.json')
ww_comb = utl.json_deserialise('/home/cbarnes/work/code/connectome/paper/real_graphs/whole_no_np/strong_only/ww'
                                    '/graph.json')

ac_mult = MultiplexConnectome(ac_comb)
ww_mult = MultiplexConnectome(ww_comb)


nodelist = sorted(ww_comb.nodes())


def filenames_gen(root_path, ext='.npy', start=0, stop=None):
    i = start - 1
    while True:
        i += 1
        if i == stop:
            raise StopIteration
        yield os.path.join(root_path, '{:03}{}'.format(i, ext))


out_root = '/home/cbarnes/work/code/connectome/paper/use_bct/controls/'


def networkx_to_und_unwe_mat(G):

    return np.array(nx.to_numpy_matrix(
        G.to_undirected(),
        nodelist,
        weight=None,
        multigraph_weight=min
    )) * (1 - np.eye(len(nodelist)))


def networkx_to_dir_unwe_mat(G):
    return np.array(nx.to_numpy_matrix(
        G,
        nodelist,
        weight=None,
        multigraph_weight=min
    )) * (1 - np.eye(len(nodelist)))


def make_controls(graphs_dict, algos_dict):
    for algo_name, algo in algos_dict.items():
        algo_root = os.path.join(out_root, algo_name)
        os.makedirs(algo_root, exist_ok=True)

        for name, G in graphs_dict.items():
            print('Processing {}'.format(name))
            graph_root = os.path.join(algo_root, name)

            os.makedirs(graph_root, exist_ok=True)

            collapsed = algo(G)
            np.save(os.path.join(algo_root, name), collapsed)
            for filename in filenames_gen(graph_root, stop=reps+1):  # +1 so that small-worldedness can be calculated
                print('Generating {}'.format(filename))
                np.save(filename, bct.randmio_und(collapsed, swap_prop)[0])


def make_semirandoms(algo_names, phys_names, extrasyn_name):
    print('Generating semi-randoms')
    for algo_name in algo_names:
        print('  ' + algo_name)
        algo_root = os.path.join(out_root, algo_name)
        for name in phys_names:
            print('    ' + name)
            shutil.copy(os.path.join(algo_root, name + '_comb.npy'), os.path.join(algo_root, name + '_comb_semi.npy'))
            real_phys_adj = np.load(os.path.join(algo_root, name + '.npy'))
            graph_root = os.path.join(algo_root, name + '_comb_semi')
            os.makedirs(graph_root, exist_ok=True)
            for randomised_ma_path, out_path in zip(
                    filenames_gen(os.path.join(algo_root, extrasyn_name), stop=reps+1),
                    filenames_gen(graph_root, stop=reps+1),
            ):
                print('      generating {}'.format(out_path))
                np.save(out_path, bct.binarize(real_phys_adj + np.load(randomised_ma_path)))


if __name__ == '__main__':

    ac = ac_mult.compose('GapJunction', 'Synapse')
    ww = ww_mult.compose('GapJunction', 'Synapse')
    ma = ac_mult['Monoamine']

    ac_comb = ac_mult.compose('GapJunction', 'Synapse', 'Monoamine')
    ww_comb = ww_mult.compose('GapJunction', 'Synapse', 'Monoamine')

    graphs = {
        'ac': ac,
        'ww': ww,
        'ma': ma,
        'ac_comb': ac_comb,
        'ww_comb': ww_comb
    }

    collapse_algos = {
        'und_bin': networkx_to_und_unwe_mat
        # 'dir_bin': networkx_to_dir_unwe_mat
    }

    make_controls(graphs, collapse_algos)

    make_semirandoms(['und_bin'], ['ac', 'ww'], 'ma')
