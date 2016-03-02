#!/bin/bash

timestamp=$(date +"%Y-%m-%d_%H:%M:%S")
zip -r "di-layer_metrics_$timestamp.zip" assortativity/ betweenness_centrality/ clustering/ degrees/ density/ \
global_efficiency/ mean_betweenness_centrality/ mean_clustering/ modularity/ path_length/ small_worldness/ \
weighted_mean_clustering/