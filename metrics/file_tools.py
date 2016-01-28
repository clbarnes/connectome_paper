import networkx as nx
import numpy as np
import connectome_utils as utl
from multiplex import MultiplexConnectome
import os
import bct
import multiprocessing as mp
try:
    from metrics.shared import set_seeds, push_exceptions
except (ImportError, SystemError):
    from shared import set_seeds, push_exceptions

set_seeds()

DATA_ROOT = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data'

REPS = 100
SWAP_PROP = 10

permutations = {
    'layers': ['GapJunction', 'Synapse', 'Monoamine', 'Neuropeptide'],
    'sources': ['ac', 'ww'],
    'ma_include_weak': [True, False],
    'directed': [False],
    'weighted': [False]
}


def networkx_to_unwe_mat(G, directed=False):
    return np.array(nx.to_numpy_matrix(
        G.to_directed() if directed else G.to_undirected(),
        sorted(G.nodes()),
        weight=None,
        multigraph_weight=min
    )) * (1 - np.eye(len(G.nodes())))


def get_original(data_root=DATA_ROOT, layers=permutations['layers'], source='ww', include_weak=False, directed=False,
        weighted=False):
    if weighted:
        raise NotImplementedError("Haven't implemented directed or weighted graphs yet")

    json_path = os.path.join(data_root, 'including_weak' if include_weak else 'strong_only', source,
                             'complete.json')

    G = utl.json_deserialise(json_path)

    C = MultiplexConnectome(G)

    if isinstance(layers, str):
        return networkx_to_unwe_mat(C[layers], directed)
    else:
        return networkx_to_unwe_mat(C.compose(*layers), directed)


def get_spec_combinations():
    outlst = []

    layer_perms = [
        ['GapJunction'],
        ['Synapse'],
        ['Monoamine'],
        ['Neuropeptide'],
        ['GapJunction', 'Synapse'],
        ['GapJunction', 'Synapse', 'Monoamine'],
        ['GapJunction', 'Synapse', 'Monoamine', 'Neuropeptide']
    ]

    for layers in layer_perms:
        for source in permutations['sources']:
            for include_weak in permutations['ma_include_weak']:
                outlst.append({
                    'layers': layers,
                    'source': source,
                    'include_weak': include_weak
                })

    return outlst


def get_root_dir_list(root):
    dirset = set()
    for spec_comb in get_spec_combinations():
        dirset.add(get_control_path(
            root, **spec_comb
        ))

    return sorted(dirset)


def setup_dirs(root, layers=permutations['layers'], source='ww', include_weak=False, directed=False,
               weighted=False):
    path = get_control_path(root, layers=permutations['layers'], source='ww', include_weak=False, directed=False,
               weighted=False)

    os.makedirs(path, exist_ok=True)


def setup_all_dirs(root):
    for dir in get_root_dir_list(root):
        os.makedirs(dir, exist_ok=True)


abbreviations = {
    'GapJunction': 'gj',
    'Synapse': 'syn',
    'Monoamine': 'ma',
    'Neuropeptide': 'np'
}


def spec_to_name(layers=permutations['layers'], source='ww', include_weak=False, directed=False,
               weighted=False):
    if directed or weighted:
        raise NotImplementedError("Haven't implemented directed or weighted graphs yet")

    path_elements = []

    layers_str = abbreviations[layers] if isinstance(layers, str) else '-'.join(sorted([abbreviations[layer] for
                                                                                        layer in layers]))
    path_elements.append(layers_str)
    if 'GapJunction' in layers or 'Synapse' in layers:
        path_elements.append(source)
    if 'Monoamine' in layers:
        path_elements.append('wk' if include_weak else 'str')

    return '_'.join(path_elements)


def get_root_path(root, layers=permutations['layers'], source='ww', include_weak=False, directed=False,
               weighted=False):
    if directed or weighted:
        raise NotImplementedError("Haven't implemented directed or weighted graphs yet")

    name = spec_to_name(layers, source, include_weak, directed, weighted)

    return os.path.join(root, name)


def get_control_path(root, layers=permutations['layers'], source='ww', include_weak=False, directed=False,
               weighted=False):
    return os.path.join(
        get_root_path(root, layers, source, include_weak, directed, weighted),
        'controls'
    )


def get_real_path(root, layers=permutations['layers'], source='ww', include_weak=False, directed=False,
               weighted=False):
    return os.path.join(
        get_root_path(root, layers, source, include_weak, directed, weighted),
        'adj.npy'
    )


def make_control(source_adj_and_filepath):
    source_adj, filepath = source_adj_and_filepath[:2]
    if len(source_adj_and_filepath) > 2:
        directed = source_adj_and_filepath[2]
    else:
        directed = False
    randomiser = bct.randmio_dir if directed else bct.randmio_und
    np.save(filepath, randomiser(source_adj, SWAP_PROP)[0])
    print('  generating {}'.format(filepath))
    return True


def filename_iter(nreps=np.inf, ext='.npy'):
    i = 0
    while i < nreps:
        yield '{:03}{}'.format(i, ext)
        i += 1


def make_controls(source_adj, out_dir, n=REPS+1, directed=False):
    with mp.Pool() as p:
        set(p.imap_unordered(
            make_control,
            ((source_adj.copy(), os.path.join(out_dir, filename), directed) for filename in filename_iter(n)),
            chunksize=int(n/mp.cpu_count())
        ))


@push_exceptions
def undi_combinations_setup():
    out_root = 'graphs'
    for comb in get_spec_combinations():
        print(comb)
        control_dir = get_control_path(out_root, **comb)
        os.makedirs(control_dir, exist_ok=True)
        real_path = get_real_path(out_root, **comb)

        adj = get_original(**comb)
        np.save(real_path, adj)

        make_controls(adj, control_dir)


@push_exceptions
def di_layers_setup():
    out_root = os.path.join('graphs', 'di_layers')

    # everything except gap junctions
    di_layers = ['Synapse', 'Neuropeptide', 'Monoamine']
    for layer_name in di_layers:
        print(layer_name)
        out_dir = os.path.join(out_root, abbreviations[layer_name])
        controls_dir = os.path.join(out_dir, 'controls')
        os.makedirs(controls_dir, exist_ok=True)
        adj = get_original(layers=layer_name, source='ac', include_weak=False, directed=True)
        np.save(os.path.join(out_dir, 'adj.npy'), adj)
        make_controls(adj, controls_dir, directed=True)

    # gap junctions
    layer_name = 'GapJunction'
    print(layer_name)
    out_dir = os.path.join(out_root, abbreviations[layer_name])
    controls_dir = os.path.join(out_dir, 'controls')
    os.makedirs(controls_dir, exist_ok=True)
    adj = get_original(layers=layer_name, source='ac', include_weak=False, directed=False)
    np.save(os.path.join(out_dir, 'adj.npy'), adj)
    make_controls(adj, controls_dir, directed=False)

    # combinations
    di_combinations_from_layers()


def collapse_list_of_arrays(lst):
    adj = lst[0]
    for layer in lst[1:]:
        adj += layer

    return bct.binarize(adj)


def di_combinations_from_layers():
    out_root = os.path.join('graphs', 'di_layers')
    for combination in [
        ('gj', 'syn'),
        ('gj', 'syn', 'ma'),
        ('gj', 'syn', 'ma', 'np')
    ]:
        print('Combining {}'.format(combination))
        out_dir = os.path.join(out_root, '-'.join(sorted(combination)))
        control_dir = os.path.join(out_dir, 'controls')
        os.makedirs(control_dir, exist_ok=True)

        real_adj = collapse_list_of_arrays([np.load(os.path.join(out_root, layer, 'adj.npy')) for layer in combination])
        np.save(os.path.join(out_dir, 'adj.npy'), real_adj)

        for filename in filename_iter(REPS+1):
            print('  generating {}'.format(filename))
            adj = collapse_list_of_arrays(
                    [np.load(os.path.join(out_root, layer, 'controls', filename)) for layer in combination]
            )
            np.save(os.path.join(control_dir, filename), adj)


if __name__ == '__main__':
    # undi_combinations_setup()
    di_layers_setup()
