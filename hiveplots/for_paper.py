from hiveplotter import HivePlot
import connectome_utils as utl
from multiplex import MultiplexConnectome
from os.path import join
# import palettable
from math import log

data_root = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data'
data_paths = dict()
weakness = 'strong_only'
source = 'ac'
data_path = join(data_root, weakness, source, 'complete.json')

G = utl.json_deserialise(data_path)
M = MultiplexConnectome(G, 'etype')

print('generating parent plot')
phys = M.compose('GapJunction', 'Synapse')
for node, degree in phys.degree().items():
    phys.node[node]['log_degree'] = log(degree) if degree else 0

for log_deg in [True, False]:
    phys_H = HivePlot(phys, config_path='config/paper.ini', order_nodes_by='log_degree' if log_deg else 'degree')
    phys_H.draw()

    # edge_colours = dict(zip(['GapJunction', 'Synapse', 'Monoamine'], palettable.colorbrewer.qualitative.Set1_3.mpl_colors))

    print('generating actual plot')
    whole = M.compose('GapJunction', 'Synapse', 'Monoamine')
    H = HivePlot(whole, parent_hiveplot=phys_H, config_path='config/paper.ini')
    H.draw()
    H.save_plot('./img/gj-syn-ma_ac-str_sort-phys{}.pdf'.format('-log' if log_deg else ''))

    print('done')
