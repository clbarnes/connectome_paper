#!/bin/bash

timestamp=$(date +"%Y-%m-%d_%H:%M:%S")
zip -r "di-layer_metric_values_$timestamp.zip" meta.txt *.csv
