from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os
import json
from seaborn import despine
plt.style.use('default')

metrics = [
    'density',
    'path_length',
    'global_efficiency',
    'mean_clustering',
    'transitivity',
    'modularity',
    'assortativity',
    'betweenness_centrality',
#    'degree_distribution',
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

    return pd.DataFrame(np.array(rows), columns=metrics)


def from_json(fname):
    with open(os.path.join(control_root, fname)) as f:
        return json.load(f)

for semirandom_comb_control in [True, False]:

    for phys_name in ['ac', 'ww']:
        fig, ax_arr = plt.subplots(3, 3)
        phys_data = from_json(phys_name + '.json')
        phys_control_data = get_metrics_from_phys_controls(phys_name)
        comb_data = from_json('{}_comb{}.json'.format(phys_name, '_semi' if semirandom_comb_control else ''))
        comb_control_data = get_metrics_from_comb_controls(phys_name, semirandom_comb_control)
        for metric_name, ax in zip(metrics, ax_arr.flatten()):
            ax.boxplot([phys_control_data[metric_name], comb_control_data[metric_name]])
            ax.scatter([1, 2], [phys_data[metric_name], comb_data[metric_name]], marker='x')
            ax.set_ylabel(metric_name)
            ax.set_xticklabels(['phys', 'comb'])
        fig.suptitle('semirandom' if semirandom_comb_control else 'random')
        plt.tight_layout()
        plt.savefig('plots/metrics/{}_{}.png'.format('semirandom' if semirandom_comb_control else 'random', phys_name),
                    dpi=150)

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

        phys_dd = np.array(from_json(phys_name + '.json')['degree_distribution'])
        comb_dd = np.array(from_json(phys_name + '_comb.json')['degree_distribution'])

        fig, ax = plt.subplots(1)

        phys_x = phys_dd[:, 0]
        phys_y = phys_dd[:, 1]
        comb_x = comb_dd[:, 0]
        comb_y = comb_dd[:, 1]

        if survival:
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

        ax.set_ylabel('survival function' if survival else 'frequency')
        ax.set_xlabel('degree')
        xlims, ylims = mk_lims([phys_x, comb_x], [phys_y, comb_y], log=True)
        ax.set_xscale('log')
        ax.set_yscale('log')
        if not survival and phys_name == 'ac':
            pass
        else:
            ax.set_xlim(xlims)
            ax.set_ylim(ylims)
        ax.legend()  # loc='lower left')
        despine(fig)

        plt.tight_layout()
        plt.savefig('plots/degdist/{}_{}_{}.png'.format('dd', phys_name, 'surv' if survival else 'freq'), dpi=150)
