#!/usr/bin/env bash

set -e

python file_tools.py
python make_metrics.py
python plot_metrics.py

zip -r plots.zip plots/di_layers
echo "Don't forget to commit the changes!"