"""
similar to syn_gj_ma_fulldeg, but combining the synaptic and gap junction layers into a single network, using green for MA and magenta for the physical layer?

Is it also possible to do one with transparency turned off?
"""
from os.path import join
import connectome_utils as utl
from multiplex import MultiplexConnectome
from hiveplotter import HivePlot
from palettable import colorbrewer as cb

CONF_PATH = 'config/phys_ma_fulldeg.ini'

DATA_ROOT = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data'
weakness = "strong_only"
source = 'ac'
data_path = join(DATA_ROOT, weakness, source, 'complete.json')

etypes = ['GapJunction', 'Synapse', 'Monoamine']

G = utl.json_deserialise(data_path)
M = MultiplexConnectome(G, 'etype')
M.sub['Physical'] = M.compose('GapJunction', 'Synapse')
for src, tgt, data in M['Physical'].edges_iter(data=True):
    data['etype'] = 'Physical'
whole = M.compose('Physical', 'Monoamine')

col_dict = {
    'Physical': (1, 0, 1),
    'Monoamine': (0, 1, 0)
}

H = HivePlot(whole,
             node_class_values=['interneuron', 'motor', 'sensory'],
             edge_category_colours=col_dict,
             edge_alpha=1,
             edge_thickness_range=(0.03, 0.03),
             config_path=CONF_PATH)
H.draw()
H.save_plot('./img/phys_ma_fulldeg-THIN.pdf'.format(source, weakness))
H.dump_config(CONF_PATH)

print('done')
