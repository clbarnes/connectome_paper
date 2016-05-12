import csv


DATA = {
    'monoamine': '/home/cbarnes/work/code/connectome/construct2/extrasyn/tgt_data/'
                 'ma_edgelist.csv',
    'monoamine_including_weak': '/home/cbarnes/work/code/connectome/construct2/extrasyn/tgt_data/'
                                'ma_edgelist_include-weak.csv',
    'neuropeptide': '/home/cbarnes/work/code/connectome/construct2/extrasyn/tgt_data/'
                    'np_edgelist.csv'
}


def csv2dict(path='node_mapping.csv', order=(1, 0)):
    with open(path) as f:
        reader = csv.reader(f)
        return {tup[order[0]]: tup[order[1]] for tup in reader}


me_to_barry = csv2dict()

headers = ['src', 'tgt', 'weight', 'transmitter', 'receptor']

for graph_name, path in DATA.items():
    edgelist = []

    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            edgelist.append((
                me_to_barry[row[0]], me_to_barry[row[1]], 1, row[2], row[3]
            ))

    edgelist = sorted(edgelist)

    with open('csvs/{}.csv'.format(graph_name.lower()), 'w') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in edgelist:
            writer.writerow(row)

