#!/usr/bin/env bash

set -e

python syn_gj_ma_physdeg.py
python syn_gj_ma_fulldeg.py
python phys_ma_fulldeg.py
python syn_gj_ma_np_physdeg.py
python syn_gj_ma_np_fulldeg.py
python ma.py

timestamp=$(date +"%Y-%m-%d_%H:%M:%S")
zip -r "hiveplots_$timestamp.zip" config/ img/