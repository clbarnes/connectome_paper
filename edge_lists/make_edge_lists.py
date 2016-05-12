import connectome_utils as utl
from multiplex import MultiplexConnectome
import csv


DATA = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data/strong_only/ac/complete.json'


def csv2dict(path='node_mapping.csv', order=(1, 0)):
    with open(path) as f:
        reader = csv.reader(f)
        return {tup[order[0]]: tup[order[1]] for tup in reader}


me_to_barry = csv2dict()

g = utl.json_deserialise(DATA)
M = MultiplexConnectome(g)

headers = ['src', 'tgt', 'weight', 'transmitter', 'receptor']

for graph_name in ['Monoamine', 'Neuropeptide']:
    graph = M[graph_name]
    edgelist = []

    for src, tgt, data in graph.edges_iter(data=True):
        edgelist.append((me_to_barry[src], me_to_barry[tgt], 1, data['transmitter'], data['receptor']))

    edgelist = sorted(edgelist)

    with open('csvs/{}.csv'.format(graph_name.lower()), 'w') as f:
        writer = csv.writer(f)
        for row in edgelist:
            writer.writerow(row)

