"""
(1) showing Syn, Gap, MA with axes ordered by physical degree
"""

from os.path import join
import connectome_utils as utl
from multiplex import MultiplexConnectome
from hiveplotter import HivePlot

data_root = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data'
data_paths = dict()
weakness = "strong_only"
source = 'ac'
data_path = join(data_root, weakness, source, 'complete.json')

G = utl.json_deserialise(data_path)

M = MultiplexConnectome(G, 'etype')

whole = M.compose('GapJunction', 'Synapse', 'Monoamine')
H = HivePlot(whole, config_path='config/syn_gj_ma_physdeg.ini')
H.draw()
H.save_plot('./img/syn_gj_ma_physdeg.pdf'.format(source, weakness))

print('done')