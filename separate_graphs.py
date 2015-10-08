from paths import data_paths, real_root
import connectome_utils as utl
from multiplex import MultiplexConnectome
import os

##################################

include_weak = 'including_weak'
source = 'ac'

G = utl.json_deserialise(data_paths[include_weak][source])
M = MultiplexConnectome(G)

utl.json_serialise(M['GapJunction'], os.path.join(real_root, 'GapJunction', source, 'graph.json'))
utl.json_serialise(M['Synapse'], os.path.join(real_root, 'Synapse', source, 'graph.json'))
utl.json_serialise(M['Monoamine'], os.path.join(real_root, 'Monoamine', include_weak, 'graph.json'))
utl.json_serialise(M['Neuropeptide'], os.path.join(real_root, 'Neuropeptide', 'graph.json'))
utl.json_serialise(M.whole, os.path.join(real_root, 'whole', include_weak, source, 'graph.json'))
utl.json_serialise(
    M.compose('GapJunction', 'Synapse', 'Monoamine'),
    os.path.join(real_root, 'whole_no_np', include_weak, source, 'graph.json')
)

#################################

include_weak = 'strong_only'
source = 'ac'

G = utl.json_deserialise(data_paths[include_weak]['ac'])
M = MultiplexConnectome(G)

utl.json_serialise(M['Monoamine'], os.path.join(real_root, 'Monoamine', include_weak, 'graph.json'))
utl.json_serialise(M.whole, os.path.join(real_root, 'whole', include_weak, source, 'graph.json'))
utl.json_serialise(
    M.compose('GapJunction', 'Synapse', 'Monoamine'),
    os.path.join(real_root, 'whole_no_np', include_weak, source, 'graph.json')
)

#################################

include_weak = 'including_weak'
source = 'ww'

G = utl.json_deserialise(data_paths[include_weak][source])
M = MultiplexConnectome(G)

utl.json_serialise(M['GapJunction'], os.path.join(real_root, 'GapJunction', source, 'graph.json'))
utl.json_serialise(M['Synapse'], os.path.join(real_root, 'Synapse', source, 'graph.json'))
utl.json_serialise(M.whole, os.path.join(real_root, 'whole', include_weak, source, 'graph.json'))
utl.json_serialise(
    M.compose('GapJunction', 'Synapse', 'Monoamine'),
    os.path.join(real_root, 'whole_no_np', include_weak, source, 'graph.json')
)

##################################

include_weak = 'strong_only'
source = 'ww'

G = utl.json_deserialise(data_paths[include_weak][source])
M = MultiplexConnectome(G)

utl.json_serialise(M.whole, os.path.join(real_root, 'whole', include_weak, source, 'graph.json'))
utl.json_serialise(
    M.compose('GapJunction', 'Synapse', 'Monoamine'),
    os.path.join(real_root, 'whole_no_np', include_weak, source, 'graph.json')
)