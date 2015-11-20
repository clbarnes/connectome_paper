from os.path import join

data_root = '/home/cbarnes/work/code/connectome/construct2/combine/tgt_data'
data_paths = dict()
for weakness in ('including_weak', 'strong_only'):
    data_paths[weakness] = dict()
    for source in ('ww', 'ac'):
        data_paths[weakness][source] = join(data_root, weakness, source, 'complete.json')

control_root = '/home/cbarnes/work/code/connectome/paper/control_graphs'
real_root = '/home/cbarnes/work/code/connectome/paper/real_graphs'