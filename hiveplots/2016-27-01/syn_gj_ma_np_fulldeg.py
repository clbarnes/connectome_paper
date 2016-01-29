"""
(4) showing Syn, Gap, MA, NP with axes ordered by whole degree
"""
from os.path import join
import connectome_utils as utl
from multiplex import MultiplexConnectome
from hiveplotter import HivePlot
from palettable import colorbrewer as cb

CONF_PATH = 'config/syn_gj_ma_np_fulldeg.ini'

DATA_ROOT = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data'
weakness = "strong_only"
source = 'ac'
data_path = join(DATA_ROOT, weakness, source, 'complete.json')

etypes = ['GapJunction', 'Synapse', 'Monoamine', 'Neuropeptide']

G = utl.json_deserialise(data_path)
M = MultiplexConnectome(G, 'etype')
whole = M.compose(*etypes)

col_dict = dict(zip(sorted(etypes), cb.qualitative.Dark2_4.mpl_colors))

H = HivePlot(whole,
             node_class_values=['interneuron', 'motor', 'sensory'],
             edge_category_colours=col_dict,
             config_path=CONF_PATH)

H.draw()
H.save_plot('./img/syn_gj_ma_np_fulldeg.pdf'.format(source, weakness))
H.dump_config(CONF_PATH)

print('done')