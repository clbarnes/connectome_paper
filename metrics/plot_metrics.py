from matplotlib import pyplot as plt
import numpy as np
import os
import json
from collections import Counter
from itertools import chain
import palettable
from scipy import stats
try:
    from metrics.file_tools import filename_iter
    from metrics.shared import push_exceptions
except (ImportError, SystemError):
    from file_tools import filename_iter
    from shared import push_exceptions

DIR = True
DATA_ROOT = '/home/cbarnes/work/code/connectome/paper/metrics/graphs{}'.format('/di_layers' if DIR else '')

TITLE_SIZE = 25
XTICKLABEL_SIZE = 25
YLABEL_SIZE = 25
XLABEL_SIZE = 25
PVALUE_SIZE = 12

# XTICKLABEL_KWARGS = {'rotation': 45, 'ha': 'right'}
XTICKLABEL_KWARGS = {}

plt.style.use('default')

SPEC_NAMES_UND = [
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

SPEC_NAMES_DIR = [
    'gj',
    'syn',
    'ma',
    'np',
    'gj-syn',
    'gj-ma-syn',
    'gj-ma-np-syn'
]

SPEC_LABELS = {
    'gj': 'GJ',
    'syn': 'Syn.',
    'ma': 'MA',
    'np': 'NP',
    'gj-syn': 'Phys.',
    'gj-ma-syn': 'Aggr.',
    'gj-ma-np-syn': 'Aggr.\nw/ NP'
}


def colour_mean(*colours):
    return np.mean(colours, axis=0)


if DIR:
    col_list = palettable.colorbrewer.qualitative.Dark2_7.mpl_colors
    # colours = dict(zip(sorted(SPEC_NAMES_DIR[:4]), col_list))
    #
    # colours['gj-syn'] = colour_mean(colours['gj'], colours['syn'])
    # colours['gj-ma-syn'] = colour_mean(colours['gj'], colours['ma'], colours['syn'])
    # colours['gj-ma-np-syn'] = colour_mean(colours['gj'], colours['ma'], colours['np'], colours['syn'])
    colours = dict(zip(sorted(SPEC_NAMES_DIR[:4]) + SPEC_NAMES_DIR[4:], col_list))
else:
    with open('colours.json') as f:
        colours = json.load(f)
    colours = {spec_name: value for spec_name, (key, value) in zip(SPEC_NAMES_UND, sorted(colours.items()))}


REPS = 100


def from_json(path, key):
    with open(path) as f:
        return json.load(f)[key]


def get_feature_for(metric_name, root):
    real = from_json(os.path.join(root, 'adj.json'), metric_name)
    control = []
    for path in (os.path.join(root, 'controls', fname) for fname in filename_iter(REPS, '.json')):
        control.append(from_json(path, metric_name))

    return real, control


def plot_metric(metric_name, plot_fn, spec_names=SPEC_NAMES_UND, filename='complete.png', directed=True,
                title=False, *args, **kwargs):
    fig, ax = plt.subplots()

    data = []

    for spec_name in spec_names:
        real, control = get_feature_for(metric_name, os.path.join(DATA_ROOT, spec_name))
        data.append((spec_name, {
            'real': real,
            'control': control
        }))

    lgd = plot_fn(ax, data, *args, **kwargs)
    # plot_fn(ax, data, *args, **kwargs)
    if title:
        ax.set_title(PLOT_TITLES[metric_name], fontsize=TITLE_SIZE)

    if plot_fn in (scatter_and_cross, bars):
        plt.ylabel(PLOT_TITLES[metric_name].lower(), fontsize=YLABEL_SIZE)
    elif plot_fn == survival:  # loglog
        plt.xlabel(PLOT_TITLES[metric_name].lower() + ' (log)', fontsize=XLABEL_SIZE)
        plt.ylabel('survival function (log)', fontsize=YLABEL_SIZE)
    elif plot_fn == survival_linlin:
        plt.xlabel(PLOT_TITLES[metric_name].lower(), fontsize=XLABEL_SIZE)
        plt.ylabel('survival function', fontsize=YLABEL_SIZE)
    elif plot_fn == survival_loglin:
        plt.xlabel(PLOT_TITLES[metric_name].lower() + ' (log)', fontsize=XLABEL_SIZE)
        plt.ylabel('survival function', fontsize=YLABEL_SIZE)
    else:
        raise ValueError('Unknown plot type for x/y labels')

    plt.tight_layout()

    plots_root = 'plots/' if not directed else 'plots/di_layers/'

    os.makedirs(plots_root + '{}'.format(metric_name), exist_ok=True)
    plt.savefig(plots_root + '{}/{}'.format(metric_name, filename), dpi=140,
                bbox_extra_artists=(lgd,) if lgd else None,
                bbox_inches='tight')
    plt.close(fig)


def get_zscores(data, absolute=False):
    reals_controls = zip([d['real'] for _, d in data], [d['control'] for _, d in data])
    if absolute:
        return [abs((real - np.mean(controls))/np.std(controls)) for real, controls in reals_controls]
    else:
        return [(real - np.mean(controls))/np.std(controls) for real, controls in reals_controls]


def data_to_pval_strs(data):
    pvals = stats.norm.sf([abs(val) for val in get_zscores(data)])*2  # 2-sided
    return ['p={:.2e}'.format(val) for val in pvals]


def box_and_cross(ax, data):
    boxes = ax.boxplot([d['control'] for _, d in data])
    ax.scatter(
            range(1, len(data)+1),
            [d['real'] for _, d in data],
            color=[colours[datum_name] for datum_name, _ in data]
    )
    scatters = [d['real'] for _, d in data]
    OFFSET_PROP = 0.05
    offset = abs(ax.get_ylim()[0] - ax.get_ylim()[1])*OFFSET_PROP
    yvals_lst = [list(boxes['whiskers'][i].get_ydata()) + list(boxes['fliers'][i].get_ydata()) + [scatters[i]]
                 for i in range(len(data))]
    txt_coords = [(i, max(yvals) + offset) for i, yvals in enumerate(yvals_lst, start=1)]
    pvals = stats.norm.sf([abs(val) for val in get_zscores(data)])*2  # 2-sided
    strs = ['p={:.2e}'.format(val) for val in pvals]
    for s, (x, y) in zip(strs, txt_coords):
        ax.text(x, y, s, horizontalalignment='center', verticalalignment='bottom', fontsize=PVALUE_SIZE)
    ax.set_xticklabels([SPEC_LABELS[spec_name] for spec_name, _ in data], fontsize=XTICKLABEL_SIZE, **XTICKLABEL_KWARGS)


def scatter_and_cross(ax, data, pvals=False):
    SCATTER_SPREAD = 0.05
    maxes = []
    xticklabels = ['']

    for x, (name, this_data) in enumerate(data, 1):
        controls = this_data['control']
        xticklabels.append(SPEC_LABELS[name])

        # xs = (np.random.random(len(controls)) - 0.5) * SCATTER_SPREAD + x
        xs = np.random.normal(loc=x, scale=SCATTER_SPREAD, size=len(controls))
        ax.scatter(xs, controls, facecolor=colours[name], lw=0, alpha=0.3, marker='o', s=25)
        ax.scatter(x, this_data['real'], facecolor=colours[name], marker='D', s=35)

        maxes.append(max(this_data['real'], max(controls)))

    if pvals:
        strs = data_to_pval_strs(data)
        OFFSET_PROP = 0.05
        offset = abs(ax.get_ylim()[0] - ax.get_ylim()[1])*OFFSET_PROP

        for x, (max_, s) in enumerate(zip(maxes, strs), 1):
            ax.text(x, max_ + offset, s, horizontalalignment='center', verticalalignment='bottom', fontsize=PVALUE_SIZE)

    ax.set_xticklabels(xticklabels, fontsize=XTICKLABEL_SIZE, **XTICKLABEL_KWARGS)


def bars(ax, data):
    width = 0.8  # must be less than 1
    barlist = ax.bar([n-width/2 for n in range(1, len(data)+1)],
                     [d['real'] for _, d in data],
                     tick_label=[spec_name for spec_name, _ in data],
                     width=width)
    for i, (spec_name, _) in enumerate(data):
        barlist[i].set_color(colours[spec_name])
    ax.set_xticks(ax.get_xticks()+width/2)
    ax.set_xticklabels([SPEC_LABELS[spec_name] for spec_name, _ in data], fontsize=XTICKLABEL_SIZE, **XTICKLABEL_KWARGS)


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


def survival_linlin(ax, data):
    return survival(ax, data, (False, False))


def survival_loglin(ax, data):
    return survival(ax, data, (False, False))



PLOT_FNS = {
    'density': bars,
    'degrees': survival,
    'path_length': scatter_and_cross,
    'global_efficiency': scatter_and_cross,
    'mean_clustering': scatter_and_cross,
    'weighted_mean_clustering': scatter_and_cross,
    'clustering': survival_linlin,
    'transitivity': scatter_and_cross,
    'modularity': scatter_and_cross,
    'assortativity': scatter_and_cross,
    'small_worldness': scatter_and_cross,
    'mean_betweenness_centrality': scatter_and_cross,
    'betweenness_centrality': survival_loglin
}

PLOT_TITLES = {
    'density': 'Density',
    'degrees': 'Degree',
    'path_length': 'Characteristic path length',
    'global_efficiency': 'Global efficiency',
    'mean_clustering': 'Mean nodewise clustering coefficient',
    'weighted_mean_clustering': 'Weighted mean nodewise clustering coefficient',
    'clustering': 'Nodewise clustering coefficient',
    'transitivity': 'Transitivity',
    'modularity': 'Maximum modularity',
    'assortativity': 'Assortativity',
    'small_worldness': 'Small world coefficient',
    'mean_betweenness_centrality': 'Mean betweenness centrality',
    'betweenness_centrality': 'Betweenness centrality'
}


def exclude_from_list(lst, *exclude_substrs):
    return [
        item for item in lst
        if not any(exclude_substr in item for exclude_substr in exclude_substrs)
    ]


@push_exceptions
def main_und():
    for metric_name, plot_fn in PLOT_FNS.items():
        print('Generating ' + metric_name)
        plot_metric(metric_name, plot_fn, SPEC_NAMES_UND, filename='complete.png')

        for not_phys_src in ['ac', 'ww']:
            print("excluding {}".format(not_phys_src))

            for not_ma_tp in ['wk', 'str']:
                print('  excluding {}'.format(not_ma_tp))
                spec_names = exclude_from_list(SPEC_NAMES_UND, not_phys_src, not_ma_tp)
                plot_metric(metric_name, plot_fn, spec_names,
                            filename='not-{}_not-{}.png'.format(not_phys_src, not_ma_tp))


@push_exceptions
def main_dir():
    for metric_name, plot_fn in PLOT_FNS.items():
        print('Generating ' + metric_name)
        plot_metric(metric_name, plot_fn, SPEC_NAMES_DIR, directed=True, filename='complete.png')

        print('  excluding np')
        spec_names = exclude_from_list(SPEC_NAMES_DIR, 'np')
        plot_metric(metric_name, plot_fn, spec_names, directed=True, filename='no-np.png')


if __name__ == "__main__":
    if DIR:
        main_dir()
    else:
        main_und()
