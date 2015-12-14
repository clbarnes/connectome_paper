import bct
from oo_attempt.file_tools import filename_iter
import numpy as np
import os
import json

data_root = 'graphs'
source = 'ac'
extrasyn_tp = 'ma'

ma_tp = 'str'
if extrasyn_tp == 'ma':
    extrasyn_tp += '_' + ma_tp

extrasyn_root = os.path.join(data_root, extrasyn_tp)

real_phys = np.load(os.path.join(data_root, 'gj-syn_{}'.format(source), 'adj.npy'))
struc, mod = bct.modularity_und(real_phys)


def gen_controls(real_phys, extrasyn_root):
    for filename in filename_iter(100):
        path = os.path.join(extrasyn_root, 'controls', filename)
        control_extrasyn = np.load(path)
        yield bct.binarize(real_phys + control_extrasyn)


def actual_random(edges):
    W = np.zeros((302, 302))
    xy = []
    while len(xy) < edges:
        pair = set(np.random.randint(0, 302, 2))
        if pair not in xy and len(pair) > 1:
            xy.append(pair)
            s = sorted(pair)
            W[s[0], s[1]] = 1
    return W + W.T


def gen_rands(real_phys, extrasyn_root):
    cont = np.load(os.path.join(extrasyn_root, 'controls', '000.npy'))
    edges = np.sum(np.triu(cont, 1))
    for _ in range(100):
        yield bct.binarize(real_phys + actual_random(edges))

real_extrasyn = np.load(os.path.join(extrasyn_root, 'adj.npy'))
real_comb = bct.binarize(real_phys + real_extrasyn)

real_comb_struc, real_comb_mod = bct.modularity_und(real_comb, kci=struc)
assert np.allclose(real_comb_struc, struc)

control_comb_mods = []
for control in gen_controls(real_phys, extrasyn_root):
    control_comb_struc, control_comb_mod = bct.modularity_und(control, kci=struc)
    assert np.allclose(control_comb_struc, struc)
    control_comb_mods.append(control_comb_mod)


random_comb_mods = []
for random in gen_rands(real_phys, extrasyn_root):
    random_comb_struc, random_comb_mod = bct.modularity_und(random, kci=struc)
    assert np.allclose(random_comb_struc, struc)
    random_comb_mods.append(random_comb_mod)

d = {
    'real_phys': mod,
    'real_comb': real_comb_mod,
    'random_comb': random_comb_mods,
    'control_comb': control_comb_mods
}

out_root = 'mod_plots/{}_{}'.format(source, extrasyn_tp)
os.makedirs(out_root, exist_ok=True)
with open(os.path.join(out_root, 'data.json'), 'w') as f:
    json.dump(d, f, indent=2, sort_keys=True)
