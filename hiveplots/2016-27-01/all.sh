#!/usr/bin/env bash

set -e

python syn_gj_ma_physdeg.py
python syn_gj_ma_fulldeg.py
python syn_gj_ma_np_physdeg.py
python ma.py
