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


for weakness in ('including_weak', 'strong_only'):
    print('Weakness: {}'.format(weakness))
    for source in ('ww', 'ac'):
        print('\tSource: {}'.format(source))
        G = utl.json_deserialise(data_paths[weakness][source])

        M = MultiplexConnectome(G, 'etype')

        whole = M.compose('GapJunction', 'Synapse', 'Monoamine')
        H = HivePlot(whole, config_path='config/whole.ini')
        H.draw()
        H.save_plot('./img/whole_{}_{}.pdf'.format(source, weakness))

print('done')