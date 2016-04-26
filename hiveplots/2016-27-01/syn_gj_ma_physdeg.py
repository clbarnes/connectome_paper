"""
(1) showing Syn, Gap, MA with axes ordered by physical degree
"""
from os.path import join
import connectome_utils as utl
from multiplex import MultiplexConnectome
from hiveplotter import HivePlot
from palettable import colorbrewer as cb

CONF_PATH = 'config/syn_gj_ma_physdeg.ini'

DATA_ROOT = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data'
weakness = "strong_only"
source = 'ac'
data_path = join(DATA_ROOT, weakness, source, 'complete.json')

etypes = ['GapJunction', 'Synapse', 'Monoamine']

G = utl.json_deserialise(data_path)
M = MultiplexConnectome(G, 'etype')
phys_degrees = M.compose(*etypes[:-1]).degree()
whole = M.compose(*etypes)
for node, data in whole.nodes_iter(data=True):
    data['physdeg'] = phys_degrees[node]


# col_dict = {
#     key: value
#     for key, value in zip(sorted(etypes+ ['Neuropeptide']), cb.qualitative.Dark2_4.mpl_colors)
#     if key in etypes
# }

col_dict = {
    'Synapse': (1, 0, 1),
    'GapJunction': (0, 0, 1),
    'Monoamine': (0, 1, 0),
}

H = HivePlot(whole,
             node_class_values=['interneuron', 'motor', 'sensory'],
             order_nodes_by='physdeg',
             edge_category_colours=col_dict,
             edge_alpha=1,
             edge_thickness_range=[0.03, 0.03],
             config_path=CONF_PATH)
H.draw()
H.save_plot('./img/syn_gj_ma_physdeg.pdf'.format(source, weakness))
H.dump_config(CONF_PATH)

print('done')
