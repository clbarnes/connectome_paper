from matplotlib import pyplot as plt
import numpy as np
import os
import json
from collections import Counter
from itertools import chain
try:
    from metrics.file_tools import filename_iter
    from metrics.shared import push_exceptions
except (ImportError, SystemError):
    from file_tools import filename_iter
    from shared import push_exceptions


plt.style.use('default')

SPEC_NAMES = [
    'gj_ac',
    'gj_ww',
    'syn_ac',
    'syn_ww',
    'ma_wk',
    'ma_str',
    'np',
    'gj-syn_ac',
    'gj-syn_ww',
    'gj-ma-syn_ac_wk',
    'gj-ma-syn_ac_str',
    'gj-ma-syn_ww_wk',
    'gj-ma-syn_ww_str',
    'gj-ma-np-syn_ac_wk',
    'gj-ma-np-syn_ac_str',
    'gj-ma-np-syn_ww_wk',
    'gj-ma-np-syn_ww_str'
]


with open('colours.json') as f:
    colours = json.load(f)
colours = {spec_name: value for spec_name, (key, value) in zip(SPEC_NAMES, sorted(colours.items()))}

REPS = 100

DATA_ROOT = '/home/cbarnes/work/code/connectome/paper/oo_attempt/graphs'


def from_json(path, key):
    with open(path) as f:
        return json.load(f)[key]


def get_feature_for(metric_name, root):
    real = from_json(os.path.join(root, 'adj.json'), metric_name)
    control = []
    for path in (os.path.join(root, 'controls', fname) for fname in filename_iter(REPS, '.json')):
        control.append(from_json(path, metric_name))

    return real, control


def plot_metric(metric_name, plot_fn, spec_names=SPEC_NAMES, filename='complete.png', *args, **kwargs):
    fig, ax = plt.subplots()

    data = []

    for spec_name in spec_names:
        real, control = get_feature_for(metric_name, os.path.join(DATA_ROOT, spec_name))
        data.append((spec_name, {
            'real': real,
            'control': control
        }))

    lgd = plot_fn(ax, data, *args, **kwargs)
    ax.set_title(metric_name)
    plt.tight_layout()

    os.makedirs('plots/{}'.format(metric_name), exist_ok=True)
    plt.savefig('plots/{}/{}'.format(metric_name, filename), dpi=140,
                bbox_extra_artists=(lgd,) if lgd else None,
                bbox_inches='tight')
    plt.close(fig)


def box_and_cross(ax, data):
    ax.boxplot([d['control'] for _, d in data])
    ax.scatter(range(1, len(data)+1), [d['real'] for _, d in data], color=[colours[metric_name] for metric_name,
                                                                                                    _ in data])
    ax.set_xticklabels([spec_name for spec_name, _ in data], rotation=90)


def bars(ax, data):
    width = 0.8  # must be less than 1
    barlist = ax.bar([n-width/2 for n in range(1, len(data)+1)],
                     [d['real'] for _, d in data],
                     tick_label=[spec_name for spec_name, _ in data],
                     width=width)
    for i, (spec_name, _) in enumerate(data):
        barlist[i].set_color(colours[spec_name])
    ax.set_xticks(ax.get_xticks()+width/2)
    ax.set_xticklabels([spec_name for spec_name, _ in data], rotation=90)


def surv(lst):
    counts = np.array(sorted(Counter(lst).items()), dtype=float)
    y = counts[:, 1]
    counts[:, 1] = 1 - np.cumsum(y)/y.sum()
    return counts


def survival(ax, data, log=(True, True)):
    for spec_name, d in data:
        real = d['real']
        control = d['control']
        surv_real = surv(real)
        ax.plot(surv_real[:, 0], surv_real[:, 1], color=colours[spec_name], label=spec_name)
        surv_cont = surv(list(chain(*control)))
        ax.plot(surv_cont[:, 0], surv_cont[:, 1],
                 color=colours[spec_name], linestyle='--', label=spec_name + ' null')

    if log[0]:
        ax.set_xscale('log')
    if log[1]:
        ax.set_yscale('log')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    return ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


PLOT_FNS = {
    'density': bars,
    'degrees': survival,
    'path_length': box_and_cross,
    'global_efficiency': box_and_cross,
    'mean_clustering': box_and_cross,
    'clustering': lambda ax, data: survival(ax, data, (False, False)),
    'transitivity': box_and_cross,
    'modularity': box_and_cross,
    'assortativity': box_and_cross,
    'mean_betweenness_centrality': box_and_cross,
    'betweenness_centrality': lambda ax, data: survival(ax, data, (True, False))
}


def exclude_from_list(lst, *exclude_substrs):
    return [
        item for item in lst
        if not any(exclude_substr in item for exclude_substr in exclude_substrs)
    ]


@push_exceptions
def main():
    for metric_name, plot_fn in PLOT_FNS.items():
        print('Generating ' + metric_name)
        plot_metric(metric_name, plot_fn, SPEC_NAMES, filename='complete.png')

        for not_phys_src in ['ac', 'ww']:
            print("excluding {}".format(not_phys_src))

            for not_ma_tp in ['wk', 'str']:
                print('  excluding {}'.format(not_ma_tp))
                spec_names = exclude_from_list(SPEC_NAMES, not_phys_src, not_ma_tp)
                plot_metric(metric_name, plot_fn, spec_names,
                            filename='not-{}_not-{}.png'.format(not_phys_src, not_ma_tp))

if __name__ == "__main__":
    main()
