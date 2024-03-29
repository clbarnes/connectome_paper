import bct
import os
import numpy as np
import json
import pandas as pd
import multiprocessing as mp
try:
    from metrics.shared import set_seeds, push_exceptions
except (ImportError, SystemError):
    from shared import set_seeds, push_exceptions

set_seeds()

SMALLWORLD_CLUSTERING = 'transitivity'  # could use

def run_algos(adj, directed=False):
    d = dict()
    d['density'] = bct.density_dir(adj)[0] if directed else bct.density_und(adj)[0]
    charpath = bct.charpath(bct.distance_bin(adj), include_infinite=False)
    d['path_length'] = charpath[0]
    d['global_efficiency'] = charpath[1]
    clustering = bct.clustering_coef_bd(adj) if directed else bct.clustering_coef_bu(adj)
    d['mean_clustering'] = clustering.mean()
    d['weighted_mean_clustering'] = np.average(clustering, weights=(adj + adj.T).sum(axis=1)),
    d['clustering'] = list(clustering)
    d['transitivity'] = bct.transitivity_bd(adj) if directed else bct.transitivity_bu(adj)  # nb in dir case,
    # transitivity is 0 a lot
    d['modularity'] = bct.modularity_dir(adj)[1] if directed else bct.modularity_und(adj)[1]
    d['assortativity'] = bct.assortativity_bin(adj + adj.T, 0)  # out-out degree correlation if directed
    betweenness_centrality = bct.betweenness_bin(adj)
    d['mean_betweenness_centrality'] = betweenness_centrality.mean()
    d['betweenness_centrality'] = list(betweenness_centrality)
    d['degrees'] = list(adj.sum(axis=1))

    return d


def dump_metrics(adj, path, directed=False):
    with open(path, 'w') as f:
        json.dump(run_algos(adj, directed), f, indent=2, sort_keys=True)


def par_dump_metrics(adj_path):
    print('generating metrics for {}'.format(adj_path))
    directed = 'di_layers' in adj_path

    adj = np.load(adj_path)
    out_path = adj_path[:-4] + '.json'
    dump_metrics(adj, out_path, directed)
    return True


def add_to_json(key, value, path):
    with open(path) as f:
        d = json.load(f)
    d[key] = value
    with open(path, 'w') as f:
        json.dump(d, f, sort_keys=True, indent=2)


def get_from_json(keys, path):
    with open(path) as f:
        d = json.load(f)
    return [d[key] for key in keys]


def get_from_jsons_in_dir(keys, root_path):
    rows = []
    for file_path in (os.path.join(root_path, filename) for filename in sorted(os.listdir(root_path))):
        if file_path.endswith('.json'):
            rows.append(get_from_json(keys, file_path))

    return pd.DataFrame(np.array(rows), columns=keys)


def mean_except(arr, ignore):
    return np.ma.array(arr, mask=[i == ignore for i in range(len(arr))]).mean()


def get_small_worldness(clustering, pathlength, rand_clustering, rand_pathlength):
    gamma = clustering/rand_clustering
    lambda_ = pathlength/rand_pathlength
    return gamma/lambda_


def get_small_worldnesses(clusterings, path_lengths):
    S = []
    for i, (clustering, path_length) in enumerate(zip(clusterings, path_lengths)):
        rand_clustering = mean_except(clusterings, i)
        rand_pathlength = mean_except(path_lengths, i)
        S.append(get_small_worldness(clustering, path_length, rand_clustering, rand_pathlength))

    return S


def add_smallworld(*spec_dirs):
    for spec_dir in spec_dirs:
        print('Generating small-worldness for {}'.format(spec_dir))
        real_path = os.path.join(spec_dir, 'adj.json')
        control_dir = os.path.join(spec_dir, 'controls')
        C_lambda = get_from_jsons_in_dir([SMALLWORLD_CLUSTERING, 'path_length'], control_dir)
        print('    controls')
        Ss = get_small_worldnesses(C_lambda[SMALLWORLD_CLUSTERING], C_lambda['path_length'])
        for val, json_name in zip(
                Ss,
                (filename for filename in sorted(os.listdir(control_dir)) if filename.endswith('.json'))
        ):
            add_to_json('small_worldness', val, os.path.join(control_dir, json_name))

        real_clustering, real_pathlength = get_from_json([SMALLWORLD_CLUSTERING, 'path_length'], real_path)
        rand_clustering = mean_except(C_lambda[SMALLWORLD_CLUSTERING], 100)
        rand_pathlength = mean_except(C_lambda['path_length'], 100)

        print('    real')
        add_to_json('small_worldness',
                    get_small_worldness(real_clustering, real_pathlength, rand_clustering, rand_pathlength),
                    real_path)


def list_adj_paths(root):
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        if 'di_layers' in dirpath and 'di_layers' not in root:
            continue
        for filename in filenames:
            path = os.path.join(dirpath, filename)
            if path.endswith('.npy'):
                out.append(path)

    return out


def make_most_metrics(root='graphs', parallel=True):
    adj_paths = list_adj_paths(root)
    if parallel:
        with mp.Pool() as p:
            set(p.imap_unordered(par_dump_metrics, adj_paths, chunksize=int(len(adj_paths) / mp.cpu_count())))
    else:
        set(map(par_dump_metrics, adj_paths))


def get_control_type_roots(root):
    return [os.path.join(root, dirn) for dirn in os.listdir(root) if os.path.isdir(os.path.join(root, dirn))]


@push_exceptions
def main():
    make_most_metrics('graphs/di_layers/')
    add_smallworld(*get_control_type_roots('graphs/di_layers/'))


if __name__ == "__main__":
    main()
