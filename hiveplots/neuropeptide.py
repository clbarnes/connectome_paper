from os.path import join
import connectome_utils as utl
from multiplex import MultiplexConnectome
from hiveplotter import HivePlot

data_root = '/home/cbarnes/code/connectome/construct2/combine/tgt_data'
data_paths = dict()
for weakness in ('including_weak', 'strong_only'):
    data_paths[weakness] = dict()
    for source in ('ww', 'ac'):
        data_paths[weakness][source] = join(data_root, weakness, source, 'complete.json')

weakness = 'strong_only'
print('Weakness: {}'.format(weakness))
source = 'ww'
print('\tSource: {}'.format(source))
G = utl.json_deserialise(data_paths[weakness][source])

M = MultiplexConnectome(G, 'etype')

np = M['Neuropeptide']
H = HivePlot(np, config_path='config/neuropeptide.ini')
H.draw()
H.save_plot('./img/neuropeptide_{}.pdf'.format(weakness))

print('done')