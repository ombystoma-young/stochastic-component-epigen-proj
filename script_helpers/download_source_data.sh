#!/bin/bash
# This script is download data for the project into a folder named source_data

mkdir source_data

# MESA methylomics peripheral CD4 (T cell)
echo "Downloading GSE56581..."
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56581/suppl/GSE56581%5Fmethylome%5Fsignal%5Fintensities%2Etxt%2Egz" -O source_data/CD4_methylome_signal_intensities.txt.gz
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56581/matrix/GSE56581_series_matrix.txt.gz" -O source_data/CD4_series_matrix.txt.gz

echo "Downloading GSE56046..."
# MESA methylomics peripheral CD14 (monocytes)
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56046/suppl/GSE56046%5Fmethylome%5Fsignal%5Fintensities%2Etxt%2Egz" -O source_data/CD14_methylome_signal_intensities.txt.gz
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56046/matrix/GSE56046_series_matrix.txt.gz" -O source_data/CD14_series_matrix.txt.gz

# Download data from Horvath, 2013
echo "Download data from Horvath's paper (S3)..."
wget -c "https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2013-14-10-r115/MediaObjects/13059_2013_3156_MOESM3_ESM.csv" -O source_data/horvath_353_cpgs_s3.csv
