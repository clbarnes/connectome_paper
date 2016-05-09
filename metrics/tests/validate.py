import numpy as np
import unittest as ut
import pandas as pd
import os
import json
from nose_parameterized import parameterized
from abc import abstractclassmethod
import csv
from datetime import datetime

DIR_CSV_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'barry_data', 'basic_measures_dir.csv')
UND_CSV_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'barry_data', 'basic_measures_und.csv')

metrics = {
    'Density': 'density',
    'Modularity': 'modularity',
    'Clustering': 'transitivity',
    'Path length': 'path_length',
    'Small worldness': 'small_worldness',
    'Assortativity': 'assortativity'
#    ('Robustness'),
#    ('Reciprocity')
}


def name_to_path(name, dir=True):
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        '..',
        'graphs',
        'di_layers' if dir else 'und_layers',
        name,
        'adj.json'
    )

graphs = [
    (' Network: Syn', 'syn'),
    (' Network: Gap', 'gj'),
    (' Network: Wired', 'gj-syn'),
    (' Network: MA', 'ma'),
    (' Network: All', 'gj-ma-syn')
]


def load_barry_data(dir=True):
    path = DIR_CSV_PATH if dir else UND_CSV_PATH
    return pd.read_csv(path, index_col=0)


def test_name(testcase_func, param_num, param):
    return '{}_{}'.format(testcase_func.__name__, param.args[1])


class BinTests:
    barry_data_path = None

    @classmethod
    def setUpClass(cls):
        cls.barry_data = load_barry_data(cls.barry_data_path)

    def assert_same_vals(self, long_metric_name, long_net_name, short_net_name):
        barry_val = self.barry_data[long_net_name][long_metric_name]
        with open(name_to_path(short_net_name)) as f:
            my_val = json.load(f)[metrics[long_metric_name]]

        if not np.allclose(barry_val, my_val, atol=1e-4):
            raise AssertionError('{} {}\n\tBarry: {}, \n\tMe: {}\n\tDiff: {}'.format(
                        dict(graphs)[long_net_name], long_metric_name.lower(), barry_val, my_val, barry_val-my_val
                    ))


    @parameterized.expand(graphs, testcase_func_name=test_name)
    def test_density(self, long_net_name, short_net_name):
        self.assert_same_vals('Density', long_net_name, short_net_name)

    @parameterized.expand(graphs, testcase_func_name=test_name)
    def test_modularity(self, long_net_name, short_net_name):
        self.assert_same_vals('Modularity', long_net_name, short_net_name)

    @parameterized.expand(graphs, testcase_func_name=test_name)
    def test_transitivity(self, long_net_name, short_net_name):
        self.assert_same_vals('Clustering', long_net_name, short_net_name)

    @parameterized.expand(graphs, testcase_func_name=test_name)
    def test_path_length(self, long_net_name, short_net_name):
        self.assert_same_vals('Path length', long_net_name, short_net_name)

    @parameterized.expand(graphs, testcase_func_name=test_name)
    def test_small_worldness(self, long_net_name, short_net_name):
        self.assert_same_vals('Small worldness', long_net_name, short_net_name)

    @parameterized.expand(graphs, testcase_func_name=test_name)
    def test_assortativity(self, long_net_name, short_net_name):
        self.assert_same_vals('Assortativity', long_net_name, short_net_name)


class BinDirTests(BinTests, ut.TestCase):
    barry_data_path = DIR_CSV_PATH


if __name__ == '__main__':
    df = load_barry_data()
