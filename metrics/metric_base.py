import networkx as nx
from abc import ABCMeta, abstractmethod, abstractproperty
import os
from multiprocessing import Pool, cpu_count
import numpy as np
from collections import defaultdict, namedtuple
import shelve
from matplotlib import pyplot as plt


GraphSpec = namedtuple('GraphSpec', ['etype', 'include_weak', 'source'])
GraphSpec.__new__.__defaults__ = (None, None, None)


class Metric():
    __metaclass__ = ABCMeta
    control_root = '/home/cbarnes/code/connectome/paper/control_graphs'
    img_root = '/home/cbarnes/code/connectome/paper/metrics/img'
    real_root = '/home/cbarnes/code/connectome/paper/real_graphs'

    @abstractproperty
    def metric_name(self):
        raise NotImplementedError

    def __init__(self, db_path, include_np=True):
        # self.multiplex = mult
        self.include_np = include_np
        self.etypes = ['GapJunction', 'Synapse', 'Monoamine']
        if include_np:
            self.etypes.append('Neuropeptide')
        self.db_path = db_path

    def generate_graph_specs(self):
        """
        Iterate through all specifications of graphs addressed, returning a namedtuple with fields for etype,
         include_weak and source (defaults None).

        :return: GraphSpec
        """
        for etype in self.etypes + ['whole']:
            if etype in ['GapJunction', 'Synapse']:
                for source in ['ac', 'ww']:
                    yield GraphSpec(etype=etype, source=source)
            if etype == 'Monoamine':
                for include_weak in ['including_weak', 'strong_only']:
                    yield GraphSpec(etype=etype, include_weak=include_weak)
            if etype == 'whole':
                for include_weak in ['including_weak', 'strong_only']:
                    for source in ['ac', 'ww']:
                        yield GraphSpec(etype=etype, include_weak=include_weak, source=source)
            if etype == 'Neuropeptide':
                yield GraphSpec(etype=etype)

    def generate_db_paths(self):
        """
        Iterate through all graph specs and return the information required to write to the database, as well as the
         graph specification

        :return: tuple((dict, str), GraphSpec)
        """
        with shelve.open(self.db_path, writeback=True) as db:
            db[self.metric_name] = tree()
            for graph_spec in self.generate_graph_specs():
                yield self.db_path_from_graph_spec(graph_spec, db), graph_spec

    def db_path_from_graph_spec(self, graph_spec, open_db):
        parent = open_db[self.metric_name]
        child_key = graph_spec.etype
        if graph_spec.include_weak:
            parent = parent[child_key]
            child_key = graph_spec.include_weak
        if graph_spec.source:
            parent = parent[child_key]
            child_key = graph_spec.source

        return parent, child_key


    @abstractmethod
    def apply_algorithm(self, graph):
        """
        Apply whatever algorithm this class is being used for to the specified graph.

        """
        raise NotImplementedError

    @abstractmethod
    def mangle_graph(self, graph):
        """
        Munge the graph into whatever form it needs to be in for the algorithm.

        """
        raise NotImplementedError

    def iterate_control_paths(self, n=100, etype=None, source=None, include_weak=None):
        """
        Iterate through the paths to the null model graphs

        :param n: Number of controls to iterate through
        :param etype: 'GapJunction', 'Synapse', 'Monoamine', 'Neuropeptide' or 'whole'
        :param source: 'ac' or 'ww'
        :param include_weak: 'including_weak' or 'strong_only'
        :return: Path to a serialised control graph
        """
        control_path = self.get_data_path(self.control_root, etype, source, include_weak)

        for filename in sorted(os.listdir(control_path))[:n]:
            yield os.path.join(control_path, filename)

    def get_data_path(self, root, etype=None, source=None, include_weak=None):
        if etype == 'whole':
            path = os.path.join(root, etype if self.include_np else etype + '_no_np', include_weak, source)
        else:
            path = os.path.join(
                root,
                *[dirname for dirname in (etype, include_weak, source) if dirname is not None]
            )

        return path

    def get_control_graph(self, path):
        raise NotImplementedError  # todo: decide on serialisation for controls

    def control_path_to_metric(self, path):
        """
        Munge and then process the graph at the location specified by path.

        """
        G = self.get_control_graph(path)
        G = self.mangle_graph(G)
        return self.apply_algorithm(G)

    def control_paths_to_metrics(self, control_paths):
        """
        Process all control graphs given in control_paths

        """
        with Pool(cpu_count()) as p:
            values = p.map(
                self.control_path_to_metric,
                control_paths,
                chunksize=int(len(control_paths)/cpu_count())
            )

        return values

    @abstractmethod
    def plot(self, etype=None, source=None, include_weak=None, graph_spec=None, show=True):
        raise NotImplementedError

    def get_data(self, etype=None, source=None, include_weak=None, graph_spec=None):
        if not graph_spec:
            graph_spec = GraphSpec(etype=etype, source=source, include_weak=include_weak)

        with shelve.open(self.db_path, flag='r') as db:
            parent, child = self.db_path_from_graph_spec(graph_spec, db)
            return parent[child]

    def plot_all(self):
        for graph_spec in self.generate_graph_specs():
            self.plot(graph_spec=graph_spec, show=False)
            plt.savefig(
                os.path.join(
                    Metric.control_root,
                    self.metric_name,
                    '_'.join(
                        *[spec for spec in (graph_spec.etype, graph_spec.include_weak, graph_spec.source)
                          if spec is not None]
                    ) + '.png'
                )
            )

    def process_graph_type(self, n_controls=100, etype=None, source=None, include_weak=None):
        """
        For a given graph type, process the real values and statistics of the null models.

        """
        if etype != 'whole':
            real_graph = self.multiplex[etype]
        else:
            real_graph = self.multiplex.compose(*self.etypes)

        real_graph = self.mangle_graph(real_graph)
        real_val = self.apply_algorithm(real_graph)

        control_paths = list(
            self.iterate_control_paths(n_controls, etype=etype, source=source, include_weak=include_weak)
        )
        control_vals = self.control_paths_to_metrics(control_paths)

        control_mean = np.mean(control_vals)
        control_sd = np.std(control_vals)

        return {'real': real_val, 'controls': control_vals, 'control_mean': control_mean, 'control_sd': control_sd}

    def process(self, n_controls=100):
        """
        For all graph types, process the real values and statistics of the null models and write to the database.

        """
        for (parent, child_key), graph_spec in self.generate_db_paths():
            parent[child_key] = self.process_graph_type(n_controls, **graph_spec.__dict__)


def tree():
    return defaultdict(tree)