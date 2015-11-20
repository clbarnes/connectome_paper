import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os
import json
from seaborn import despine
from collections import Counter
from itertools import chain
plt.style.use('default')

metrics = [
    'density',
    'path_length',
    'global_efficiency',
    'clustering',
    'transitivity',
    'modularity',
    'assortativity',
    'betweenness_centrality',
    # 'degrees',
    'small_worldness'
]

control_root = '/home/cbarnes/work/code/connectome/paper/use_bct/controls/und_bin/'


def filenames_gen(root_path, ext='.npy', start=0, stop=None):
    i = start - 1
    while True:
        i += 1
        if i == stop:
            raise StopIteration
        yield os.path.join(root_path, '{:03}{}'.format(i, ext))


def get_metrics_from_phys_controls(name):
    rows = []
    for path in filenames_gen(os.path.join(control_root, name), ext='.json', stop=100):
        with open(path) as f:
            d = json.load(f)
        rows.append([d[metric] for metric in metrics])

    return pd.DataFrame(np.array(rows), columns=metrics)


def get_metrics_from_comb_controls(name, semirandom=True):
    rows = []
    for path in filenames_gen(
            os.path.join(control_root, '{}_comb{}'.format(name, '_semi' if semirandom else '')),
            ext='.json', stop=100
    ):
        with open(path) as f:
            d = json.load(f)
        rows.append([d[metric] for metric in metrics])

    return pd.DataFrame(np.array(rows), columns=metrics, dtype='object')


def from_json(fname):
    with open(os.path.join(control_root, fname)) as f:
        return json.load(f)


def box_and_cross(ax: matplotlib.axes._subplots.AxesSubplot, phys_real, comb_real, phys_control_column,
                  comb_control_column):
    ax.boxplot([phys_control_column, comb_control_column])
    ax.scatter([1, 2], [phys_real, comb_real], marker='x')
    ax.set_xticklabels(['phys', 'comb'])
    # ax.set_ylabel(metric_name)


def bars(ax: matplotlib.axes._subplots.AxesSubplot, phys_real, comb_real):
    width = 0.8
    ax.bar([n+1-0.5*width for n in [1, 2]], [phys_real, comb_real], width=width)
    ax.set_xticklabels(['phys', 'comb'])


def points2dist(points):
    """
    Value=x, freq=y
    """

    counts = Counter(points)
    xy_pairs = []
    for x, y in sorted(counts.items()):
        xy_pairs.append([x, y])

    return np.array(xy_pairs)


def deg_dist(ax: matplotlib.axes._subplots.AxesSubplot, phys_real, comb_real, plot_survival=True):
        phys_dd = points2dist(phys_real)
        comb_dd = points2dist(comb_real)

        phys_x = phys_dd[:, 0]
        phys_y = phys_dd[:, 1]
        comb_x = comb_dd[:, 0]
        comb_y = comb_dd[:, 1]

        if plot_survival:
            phys_y, comb_y = surv(phys_y), surv(comb_y)
            ax.plot(phys_x, phys_y, color='b', label='phys')
            ax.plot(comb_x, comb_y, color='r', label='comb')

            phys_e = np.average(phys_x, weights=phys_y)
            ax.plot(*bestfit_log(phys_x[phys_x > phys_e], phys_y[phys_x > phys_e]), color='k', linestyle='--')
            comb_e = np.average(comb_x, weights=comb_y)
            ax.plot(*bestfit_log(comb_x[comb_x > comb_e], comb_y[comb_x > comb_e]), color='k', linestyle='--')
        else:
            ax.scatter(phys_x, phys_y, marker='x', color='b', label='phys')
            ax.scatter(comb_x, comb_y, marker='x', color='r', label='comb')

        ax.set_ylabel('survival function' if plot_survival else 'frequency')
        ax.set_xlabel('degree')
        ax.set_xscale('log')
        ax.set_yscale('log')
        # xlims, ylims = mk_lims([phys_x, comb_x], [phys_y, comb_y], log=True)
        # if not survival and phys_name == 'ac':
        #     pass
        # else:
        #     ax.set_xlim(xlims)
        #     ax.set_ylim(ylims)
        ax.legend()


def surv_plots(ax: matplotlib.axes._subplots.AxesSubplot, phys_real, phys_control_column, comb_real, comb_control_column):
    phys_dist = points2dist(phys_real)
    comb_dist = points2dist(comb_real)
    phys_control_dist = points2dist(list(chain(phys_control_column)))
    comb_control_dist = points2dist(list(chain(comb_control_column)))

    ax.plot(phys_dist[:, 0], surv(phys_dist[:, 1]), color='b', label='phys')
    ax.plot(phys_control_dist[:, 0], surv(phys_control_dist[:, 1]), color='b', linestyle='--', label='phys rand')
    ax.plot(comb_dist[:, 0], surv(comb_dist[:, 1]), color='r', label='comb')
    ax.plot(comb_control_dist[:, 0], surv(comb_control_dist[:, 1]), color='r', linestyle='--', label='comb rand')

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.legend()


global_plots = [
    'density',
    'path_length',
    'global_efficiency',
    # 'clustering',
    'mean_clustering',
    'transitivity',
    'modularity',
    'assortativity',
    # 'betweenness_centrality',
    'mean_betweenness_centrality',
    # 'degrees',
    'small_worldness'
]


## GLOBAL PLOTS

for semirandom_comb_control in [True, False]:

    for phys_name in ['ac', 'ww']:
        fig, ax_arr = plt.subplots(3, 3)
        phys_data = from_json(phys_name + '.json')
        phys_control_data = get_metrics_from_phys_controls(phys_name)
        comb_data = from_json('{}_comb{}.json'.format(phys_name, '_semi' if semirandom_comb_control else ''))
        comb_control_data = get_metrics_from_comb_controls(phys_name, semirandom_comb_control)
        for metric_name, ax in zip(global_plots, ax_arr.flatten()):
            box_and_cross(ax, phys_data, comb_data, phys_control_data, comb_control_data)
            ax.set_ylabel(metric_name)
        fig.suptitle('semirandom' if semirandom_comb_control else 'random')
        plt.tight_layout()
        plt.savefig('plots/global_metrics/{}_{}.png'.format('semirandom' if semirandom_comb_control else 'random',
                                                       phys_name), dpi=150)

# degree_distribution


def surv(arr):
    return 1 - np.cumsum(arr)/arr.sum()


def mk_lim(datas, prop, log=False):
    xmin = np.min([np.min(x) for x in datas])
    xmax = np.max([np.max(x) for x in datas])
    range_ = xmax - xmin
    llim, ulim = xmin-range_*prop, xmax + range_*prop
    return llim if not log else max(llim, xmin.min()), ulim


def mk_lims(x_datas, y_datas, prop=0.3, log=False):
    return mk_lim(x_datas, prop, log), mk_lim(y_datas, prop, log)


def bestfit(x, y, deg=1):
    fit = np.polyfit(x, y, deg)
    return x, np.polyval(fit, x)


def bestfit_log(x, y, deg=1):
    nonzeros = x*y > 0
    fit = np.polyfit(np.log(x[nonzeros]), np.log(y)[nonzeros], deg)
    print(fit)
    return x[nonzeros], np.exp(np.polyval(fit, np.log(x[nonzeros])))


for survival in [True, False]:
    for phys_name in ['ac', 'ww']:

        phys_degrees = np.array(from_json(phys_name + '.json')['degrees'])
        comb_degrees = np.array(from_json(phys_name + '_comb.json')['degrees'])

        fig, ax = plt.subplots(1)

        deg_dist(ax, phys_degrees, comb_degrees, plot_survival=survival)

        despine(fig)

        plt.tight_layout()
        plt.savefig('plots/nodewise_metrics/{}_{}_{}.png'.format('degdist', phys_name, 'surv' if survival else 'freq'),
                    dpi=150)
