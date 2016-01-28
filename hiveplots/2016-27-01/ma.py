"""
(4) showing the various monoamines only (with a different set of colours to the other plots)
"""
from os.path import join
import connectome_utils as utl
from multiplex import MultiplexConnectome
from hiveplotter import HivePlot
from palettable import colorbrewer as cb

DATA_ROOT = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data'
CONF_PATH = 'config/ma.ini'
weakness = "strong_only"
source = 'ac'
data_path = join(DATA_ROOT, weakness, source, 'complete.json')

G = utl.json_deserialise(data_path)
M = MultiplexConnectome(G, 'etype')
whole = M['Monoamine']

transmitters = set(data['transmitter'] for src, tgt, data in whole.edges_iter(data=True))

col_dict = dict(zip(sorted(transmitters), cb.qualitative.Set1_4.mpl_colors))

H = HivePlot(whole,
             edge_category_colours=col_dict,
             config_path=CONF_PATH)
H.draw()
H.save_plot('./img/ma.pdf'.format(source, weakness))
H.dump_config(CONF_PATH)

print('done')