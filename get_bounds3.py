import networkx as nx
import connectome_utils as utl
from multiplex import MultiplexConnectome
from collections import namedtuple
import csv
import logging
import progressbar as pb

"""
Generate some statistics about the extrasynaptic connectomes.

Monoamine:
    - In total
    - By monoamine
    - By receptor?

Neuropeptide:
    - In total
    - By neuropeptide
    - By receptor

Extrasynaptic broadcaster
- how many nodes express a monoamine, but do not synapse to a node with a receptor for that monoamine
Extrasynaptic receiver
- how many nodes express a monoamine receptor, but do not receive a synapse from a node expressing that monoamine
Possibly aminergic synapses from broadcasters
- how many synapses made by monoamine-expressing nodes are to cognate receptor-expressing nodes
Possibly aminergic synapses
- how many synapses may be aminergic

"""

FORMAT = '%(asctime) -- %(message)'
logging.basicConfig(format=FORMAT, filename='logs.txt', filemode='w')

# sources = ['ac', 'ww']
sources = ['ac']
etypes = {'ma': 'Monoamine', 'np': 'Neuropeptide'}

DATA_SRC = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data/strong_only/ac/complete.json'

TOTAL_NODES = 302

Edge = namedtuple('ExtrasynapticEdge', ['src', 'tgt', 'transmitter', 'receptor'])


def get_transmitters(etype):
    G = utl.json_deserialise(DATA_SRC)
    out = set()
    for src, tgt, data in G.edges_iter(data=True):
        if data['etype'] == etypes[etype]:
            out.add(data['transmitter'])

    return sorted(out)


def get_receptors(etype):
    G = utl.json_deserialise(DATA_SRC)
    out = set()
    for src, tgt, data in G.edges_iter(data=True):
        if data['etype'] == etypes[etype]:
            out.add(data['receptor'])

    return sorted(out)


class NodeInfo():
    def __init__(self, etype='ma', source='ac', transmitter=None, receptor=None):
        # assert None in [transmitter, receptor], 'One of transmitter and receptor must be None'
        self.source = source
        self.etype = etype
        self.transmitter = transmitter
        self.receptor = receptor
        self.extrasyn = self.filtered_extrasyn()
        self.syn = self.get_syn_graph()
        self._info = self._make_info()

    def get_syn_graph(self):
        return MultiplexConnectome(utl.json_deserialise(DATA_SRC))['Synapse']

    def filtered_extrasyn(self):
        G = utl.json_deserialise(DATA_SRC)
        extrasyn = MultiplexConnectome(G)[etypes[self.etype]]
        filtered = nx.MultiDiGraph()
        filtered.add_nodes_from(G.nodes())

        for src, tgt, data in extrasyn.edges_iter(data=True):
            if self.transmitter and self.transmitter != data['transmitter']:
                continue
            if self.receptor and self.receptor != data['receptor']:
                continue

            filtered.add_edge(src, tgt, attr_dict=data)

        return filtered

    def _make_info(self):
        node_info = dict()

        for src, tgt, data in self.extrasyn.edges_iter(data=True):
            transmitter = data['transmitter']
            receptor = data['receptor']

            for node in (src, tgt):
                if node not in node_info:
                    node_info[node] = {'sends': set(), 'receives': set(), 'receptors': set()}

            node_info[src]['sends'].add(transmitter)
            node_info[tgt]['receives'].add(transmitter)
            node_info[tgt]['receptors'].add(receptor)

        return node_info

    def __getitem__(self, item):
        return self._info.get(item, {'sends': set(), 'receives': set(), 'receptors': set()})

    def __iter__(self):
        return self._info.__iter__()

    def items(self):
        return self._info.items()

    def keys(self):
        return self._info.keys()

    def values(self):
        return self._info.values()

    def __len__(self):
        return len(self._info)


def to_pc(numerator, denominator, raw=False):
    if raw:
        return '{} of {}, {:.2f}%'.format(numerator, denominator, numerator/denominator * 100)
    else:
        return numerator/denominator


class NodeInfoStatistics():
    def __init__(self, info):
        self.info = info

    def get_extrasynaptic_broadcasters(self, raw=False):
        """
        Proportion of all nodes expressing the transmitter, which do not synapse to nodes expressing a cognate receptor.

        :return: float
        """

        broadcasters = set()
        extrasyn_broadcasters = set()

        for extrasyn_src, data in self.info.items():
            if not data['sends']:
                continue

            if self.info.transmitter is None or self.info.transmitter in data['sends']:
                broadcasters.add(extrasyn_src)
            else:
                continue

            for transmitter in data['sends']:
                is_extrasyn_broadcaster = True
                for syn_tgt in self.info.syn.edge[extrasyn_src]:
                    if transmitter in self.info[syn_tgt]['receives']:
                        is_extrasyn_broadcaster = False
                        break
                if is_extrasyn_broadcaster:
                    extrasyn_broadcasters.add(extrasyn_src)

        return to_pc(len(extrasyn_broadcasters), len(broadcasters), raw)

    def get_extrasynaptic_receivers(self, raw=False):
        """
        Proportion of all nodes expressing the receptor, which do not synapse to nodes expressing the cognate transmitter.

        :return: float
        """

        receivers = set()
        extrasyn_receivers = set()

        syn_rev = nx.reverse(self.info.syn)
        for extrasyn_tgt, data in self.info.items():
            if not data['receives']:
                continue

            if self.info.receptor is None or self.info.receptor in data['receptors']:
                receivers.add(extrasyn_tgt)
            else:
                continue

            for transmitter in data['receives']:
                is_extrasyn_receiver = True
                for syn_src in syn_rev.edge[extrasyn_tgt]:
                    if transmitter in self.info[syn_src]['sends']:
                        is_extrasyn_receiver = False
                        break
                if is_extrasyn_receiver:
                    extrasyn_receivers.add(extrasyn_tgt)

        return to_pc(len(extrasyn_receivers), len(receivers), raw)

    def get_extrasynaptic_in_parallel_with_synapse(self, raw=False):
        """
        Proportion of 'extrasynaptic' edges which are actually in parallel with synapses.

        :param percentage:
        :return:
        """
        possible_synaptic = []

        syn_edges = set(self.info.syn.edges_iter())
        for src, tgt in self.info.extrasyn.edges_iter():
            if (src, tgt) in syn_edges:
                possible_synaptic.append((src, tgt))

        return to_pc(len(possible_synaptic), self.info.extrasyn.number_of_edges(), raw)

    def get_synapse_from_broadcaster_in_parallel_with_extrasyn(self, raw=False):
        """
        Proportion of synaptic edges from broadcasters which are in parallel with monoamine edges

        :param percentage:
        :return:
        """

        extrasyn_edgeset = set(self.info.extrasyn.edges_iter())

        syn_edges_from_broadcasters = []
        in_parallel = []

        for src, tgt, data in self.info.syn.edges_iter(data=True):
            if self.info.extrasyn.edge[src]:
                syn_edges_from_broadcasters.append((src, tgt))
                if (src, tgt) in extrasyn_edgeset:
                    in_parallel.append((src, tgt))

        return to_pc(len(in_parallel), len(syn_edges_from_broadcasters), raw)


def stats_to_list(row_title: str, stats_obj: NodeInfoStatistics, raw: bool=False):
    """

    :param row_title:
    :param stats_obj: NodeInfoStatistics
    :return: list
    """
    return [
        row_title,
        stats_obj.get_extrasynaptic_broadcasters(raw),
        stats_obj.get_extrasynaptic_receivers(raw),
        stats_obj.get_extrasynaptic_in_parallel_with_synapse(raw),
        stats_obj.get_synapse_from_broadcaster_in_parallel_with_extrasyn(raw)
    ]


class NullProgressBar():
    def __init__(self):
        self.value = 0

    def update(self, i):
        pass


def make_table(source, etype, raw=False, pbar=NullProgressBar()):
    with open('output_table_{}_{}_{}.tsv'.format(source, etype, 'raw' if raw else 'percent'), 'w') as f:
        w = csv.writer(f, delimiter='\t')
        headers = [
            source.upper(),
            'Extrasynaptic Broadcasters',
            'Extrasynaptic Receivers',
            'Extrasynaptic edges in parallel with synapses',
            'Synaptic edges from broadcasters which are to receivers'
        ]
        w.writerow(headers)

        # whole
        stats = NodeInfoStatistics(NodeInfo(etype=etype, source=source))
        row = stats_to_list('Whole', stats, raw)
        w.writerow(row)
        pbar.update(pbar.value + 1)

        # transmitters
        w.writerow(['TRANSMITTERS'])
        for transmitter in get_transmitters(etype):
            stats = NodeInfoStatistics(NodeInfo(etype=etype, source=source, transmitter=transmitter))
            row = stats_to_list(transmitter, stats, raw)
            w.writerow(row)
            pbar.update(pbar.value + 1)

        # receptors
        w.writerow(['RECEPTORS'])
        for receptor in get_receptors(etype):
            stats = NodeInfoStatistics(NodeInfo(etype=etype, source=source, receptor=receptor))
            row = stats_to_list(receptor, stats, raw)
            w.writerow(row)
            pbar.update(pbar.value + 1)


def count_iterations():
    to_sum = []
    for etype in etypes:
        to_sum.append(len(get_receptors(etype)))
        to_sum.append(len(get_transmitters(etype)))
        to_sum.append(1)

    return sum(to_sum)


def make_tables():
    i_per_source_per_raw = count_iterations()

    iterations = i_per_source_per_raw * len(sources) * 2

    pbar = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(), pb.ETA()], maxval=iterations).start()

    for source in sources:
        for etype in etypes:
            for raw in (True, False):
                make_table(source, etype, raw, pbar)


def main():
    make_tables()

if __name__ == '__main__':
    main()
    logging.shutdown()