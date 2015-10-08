from randomise_algorithms import *
import connectome_utils as utl
import multiprocessing as mp
import os
from paths import control_root, real_root


class ControlFactory():
    def __init__(self, graph, out_dir, algorithm):
        self.graph = graph
        self.out_dir = out_dir
        self.algorithm = algorithm

    def run(self, nreps=100):
        cores = mp.cpu_count()
        with mp.Pool(cores) as p:
            p.map(
                algorithm_wrapper,
                [(self.algorithm, self.graph.copy(), '{}/{:03d}.json'.format(self.out_dir, i)) for i in range(nreps)],
                chunksize=int(nreps/cores)
            )

short_run = True

control_dir_paths = []
real_graph_paths = []
for root, dirnames, filenames in os.walk(real_root):
    for filename in filenames:
        real_graph_paths.append(os.path.join(root, filename))
        control_dir_paths.append(root.replace(real_root, control_root))

nreps = 100

for i, (graph_path, out_dir) in enumerate(zip(real_graph_paths, control_dir_paths)):
    print('Processing:\n\t{}\n\t{} of {} graph specs'.format(
        graph_path.replace(real_root, ''),
        i+1,
        len(control_dir_paths)
    ))

    if short_run and len(os.listdir(out_dir)) == nreps:
        print('Files already exist - SKIPPING')
        continue
    ControlFactory(utl.json_deserialise(graph_path), out_dir, di_maslov).run(nreps)