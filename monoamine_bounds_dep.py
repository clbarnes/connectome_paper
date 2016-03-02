import connectome_utils as utl
from multiplex import MultiplexConnectome
from itertools import product
import pandas as pd

DATA_SRC = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data/strong_only/ac/complete.json'

G = utl.json_deserialise(DATA_SRC)


def edges_of_type(graph, *types):
    for src, tgt, data in graph.edges_iter(data=True):
        if data['etype'] in types:
            yield src, tgt, data

def edges_between(G, srcs, tgts):
    for src, tgt in product(srcs, tgts):
        if src in G.edge and tgt in G.edge[src]:
            yield src, tgt, G.edge[src][tgt]


senders = {}
receivers = {}

for src, tgt, data in edges_of_type(G, 'Monoamine'):
    try:
        senders[data['transmitter']].add(src)
    except KeyError:
        senders[data['transmitter']] = {src}

    try:
        receivers[data['transmitter']].add(tgt)
    except KeyError:
        receivers[data['transmitter']] = {tgt}

    try:
        G.node[src]['sends'].add(data['transmitter'])
    except KeyError:
        G.node[src]['sends'] = {data['transmitter']}

    try:
        G.node[tgt]['receives'].add(data['transmitter'])
    except KeyError:
        G.node[tgt]['receives'] = {data['transmitter']}


extrasyn_senders = {key: set() for key in senders}
extrasyn_receivers = {key: set() for key in receivers}

M = MultiplexConnectome(G)

for transmitter in senders:
    for src in senders[transmitter]:
        if not list(edges_between(M['Synapse'], [src], receivers[transmitter])):
            extrasyn_senders[transmitter].add(src)

    for tgt in receivers[transmitter]:
        if not list(edges_between(M['Synapse'], senders[transmitter], [tgt])):
            extrasyn_receivers[transmitter].add(tgt)

headers = ['broacasters', 'extrasyn_broadcasters', 'extrasyn_broadcast_prop', 'receivers', 'extrasyn_receivers',
           'extrasyn_receive_prop']
data = []
index = []

for transmitter in senders:
    index.append(transmitter)
    broadcaster_n = len(senders[transmitter])
    extrasyn_broadcasters = len(extrasyn_senders[transmitter])
    extrasyn_broadcast_prop = extrasyn_broadcasters / broadcaster_n * 100
    receiver_n = len(receivers[transmitter])
    extrasyn_receiver_n = len(extrasyn_receivers[transmitter])
    extrasyn_receive_prop = extrasyn_receiver_n / receiver_n * 100

    data.append([broadcaster_n, extrasyn_broadcasters, extrasyn_broadcast_prop, receiver_n, extrasyn_receiver_n,
                 extrasyn_receive_prop])

df = pd.DataFrame(data, index=index, columns=headers)


print('done')
