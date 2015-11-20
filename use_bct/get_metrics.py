import connectome_utils as utl
import bct
import os
import numpy as np
from collections import Counter
import json
import pandas as pd

reps = 100
swap_prop = 10

ww_comb = utl.json_deserialise('/home/cbarnes/work/code/connectome/paper/real_graphs/whole_no_np/strong_only/ww'
                                    '/graph.json')

nodelist = sorted(ww_comb.nodes())


def filenames_gen(root_path, ext='.npy', start=0):
    while True:
        yield os.path.join(root_path, '{:03}{}'.format(start, ext))
        start += 1


control_root = '/home/cbarnes/work/code/connectome/paper/use_bct/controls/'


def control_dir(graph_name):
    return os.path.join(control_root, graph_name)


def degree_dist(adj):
    """
    There are y nodes with a degree of x
    """

    counts = Counter(np.array(adj).sum(axis=1))
    xy_pairs = []
    for x, y in sorted(counts.items()):
        xy_pairs.append([x, y])

    return xy_pairs


def run_algos(adj):
    d = dict()
    d['density'] = bct.density_und(adj)[0]
    charpath = bct.charpath(bct.distance_bin(adj))
    d['path_length'] = charpath[0]
    d['global_efficiency'] = charpath[1]
    clustering = bct.clustering_coef_bu(adj)
    d['mean_clustering'] = clustering.mean()
    d['clustering'] = list(clustering)
    d['transitivity'] = bct.transitivity_bu(adj)
    d['modularity'] = bct.modularity_und(adj)[1]
    d['assortativity'] = bct.assortativity_bin(adj, 0)
    betweenness_centrality = bct.betweenness_bin(adj)
    d['mean_betweenness_centrality'] = betweenness_centrality.mean()
    d['betweenness_centrality'] = list(betweenness_centrality)
    # d['degree_distribution'] = degree_dist(adj)
    d['degrees'] = list(adj.sum(axis=1))

    return d


def dump_metrics(adj, path):
    with open(path, 'w') as f:
        json.dump(run_algos(adj), f, indent=2, sort_keys=True)


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


def generate_most_metrics(control_root, ignore_if_exists=False):
    """
    All except small worldness
    :param control_root:
    :return:
    """
    print('Generating metrics...')
    for dirpath, _, filenames in os.walk(control_root):
        for filename in sorted(filenames):
            if filename.endswith('.npy'):
                source_file = os.path.join(dirpath, filename)
                target_file = source_file[:-4] + '.json'
                print('  Processing {}'.format(target_file))
                if os.path.exists(target_file) and ignore_if_exists:
                    continue
                else:
                    dump_metrics(np.load(source_file), target_file)


def add_smallworld(*control_type_roots):
    for inner_control_root in control_type_roots:
        print('Generating small-worldness for {}'.format(inner_control_root))
        dir_contents = [os.path.join(inner_control_root, obj_name) for obj_name in os.listdir(inner_control_root)]
        for dirpath in (obj for obj in sorted(dir_contents) if os.path.isdir(obj)):
            print('  for {}'.format(dirpath))
            real_file = dirpath + '.json'
            C_lambda = get_from_jsons_in_dir(['mean_clustering', 'path_length'], dirpath)
            print('    controls')
            Ss = get_small_worldnesses(C_lambda['mean_clustering'], C_lambda['path_length'])
            for val, json_name in zip(
                    Ss,
                    (filename for filename in sorted(os.listdir(dirpath)) if filename.endswith('.json'))
            ):
                add_to_json('small_worldness', val, os.path.join(dirpath, json_name))

            real_clustering, real_pathlength = get_from_json(['mean_clustering', 'path_length'], real_file)
            rand_clustering = mean_except(C_lambda['mean_clustering'], 100)
            rand_pathlength = mean_except(C_lambda['path_length'], 100)

            print('    real')
            add_to_json('small_worldness',
                        get_small_worldness(real_clustering, real_pathlength, rand_clustering, rand_pathlength),
                        real_file)


if __name__ == '__main__':
    types = ['und_bin']
    generate_most_metrics(control_root)
    add_smallworld(*[os.path.join(control_root, tp) for tp in types])
