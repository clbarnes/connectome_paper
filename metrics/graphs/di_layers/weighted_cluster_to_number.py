import os
import json

def var_replace(path):
    with open(path) as f:
        data = json.load(f)
    # data['weighted_mean_clustering'] = data['weighted_mean_clustering'][0]
    with open(path, 'w') as f:
        json.dump(data, f, sort_keys=True, indent=2)

for root, dirnames, filenames in os.walk('.'):
    for filename in filenames:
        if filename.endswith('.json'):
            var_replace(os.path.join(root, filename))

