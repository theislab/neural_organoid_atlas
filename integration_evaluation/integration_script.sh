#!/bin/bash

# Has to be executed before scanvi
python scvi_21_05.py

levels=("1" "12" "123")

# Loop over the annotation levels
for level in "${levels[@]}"
do
    python scanvi_21_05.py --cell_type_key=snapseed_pca_rss_level_${level} --output_filename=scanvi_results_level_${level}
    python scpoli_21_05.py --cell_type_keys=snapseed_pca_rss_level_${level} --output_filename=scpoli_results_level_${level}
done
