from os.path import join

data_root = '/home/cbarnes/code/connectome/construct2/combine/tgt_data'
data_paths = dict()
for weakness in ('include_weak', 'strong_only'):
    data_paths[weakness] = dict()
    for source in ('ww', 'ac'):
        data_paths[weakness][source] = join(data_root, weakness, source)