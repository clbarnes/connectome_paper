from os.path import join
import connectome_utils as utl
from multiplex import MultiplexConnectome
from hiveplotter import HivePlot

data_root = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data'
data_paths = dict()
for weakness in ('including_weak', 'strong_only'):
    data_paths[weakness] = dict()
    for source in ('ww', 'ac'):
        data_paths[weakness][source] = join(data_root, weakness, source, 'complete.json')


for weakness in ('including_weak', 'strong_only'):
    print('Weakness: {}'.format(weakness))
    source = 'ww'
    print('\tSource: {}'.format(source))
    G = utl.json_deserialise(data_paths[weakness][source])

    M = MultiplexConnectome(G, 'etype')

    ma = M['Monoamine']
    H = HivePlot(ma, config_path='config/monoamine.ini')
    H.draw()
    H.save_plot('./img/monoamine_{}.pdf'.format(weakness))

print('done')