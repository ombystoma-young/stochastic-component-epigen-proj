#!/bin/bash
# This script runs alls steps of analysis

# create and activate environment
mamba env create -f env.yaml
mamba activate cba_fin_proj

# download source data
sh script_helpers/download_source_data.sh

# process source data
mkdir data

# copy Horvath's clocks weights
cp source_data/horvath_353_cpgs_s3.csv data/horvath_353_cpgs_s3.csv

# extract Horvath's CpG ids
python script_helpers/parse_cpg_positions.py -i source_data/horvath_353_cpgs_s3.csv -o data/horvaths.pkl

# find methylation fractions
python script_helpers/process_source_data.py -i source_data/CD14_methylome_signal_intensities.txt.gz \
-d source_data/CD14_series_matrix.txt.gz \
-o data/CD14_DNAm_frac.pkl -od data/CD14_sample.pkl \
-f data/horvaths.pkl
python script_helpers/process_source_data.py -i source_data/CD4_methylome_signal_intensities.txt.gz \
-d source_data/CD4_series_matrix.txt.gz \
-o data/CD4_DNAm_frac.pkl -od data/CD4_sample.pkl \
-f data/horvaths.pkl

# find optimized gamma and sigma
python3 script_helpers/search_optimal_parameters.py

# simulate artificial cohorts
python3 script_helpers/simulate_artificial_cohorts.py