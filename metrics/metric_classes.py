from .metric_base import Metric, GraphSpec
import networkx as nx
from collections import Counter
import numpy as np
from matplotlib import pyplot as plt
import shelve


class DegreeDistribution(Metric):
    def __init__(self, mult, db_path, include_np=True):
        super().__init__(mult, db_path, include_np)

    @property
    def metric_name(self):
        return 'DegreeDistribution'

    def mangle_graph(self, graph):
        return nx.Graph(graph)

    def apply_algorithm(self, graph):
        degrees = graph.degree().values()
        counted = Counter(degrees)
        return np.array(sorted([(value, key) for key, value in counted.items()], key=lambda x: x[0]))

    def plot(self, etype=None, source=None, include_weak=None, graph_spec=None, show=True):
        data = self.get_data(etype, source, include_weak, graph_spec)
        plt.figure()
        ax = plt.gca()
        ax.scatter(data)
        ax.set_xscale('log')
        ax.set_yscale('log')
        if show:
            plt.show()
