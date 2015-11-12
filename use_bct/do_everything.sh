#!/usr/bin/env bash

set -e
python mk_controls.py
python get_metrics.py
# python mk_plots.py