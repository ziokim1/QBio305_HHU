#!/bin/bash
# Change permissions on your computer so that you can run a shell script by typing: 'chmod +x setup.sh' (without the quotes) at the terminal prompt 
# Then type './setup.sh' (without the quotes) at the prompt.

# setup working directories

mkdir -p ~/Documents/QBIO305_group1_assignment_data

# download all data used for analysis
wget -nc -P ~/Documents/QBIO305_group1_assignment_data https://raw.githubusercontent.com/ziokim1/QBio305_HHU/main/data/305_ANOVA_QTL_analysis_data.csv
wget -nc -P ~/Documents/QBIO305_group1_assignment_data https://raw.githubusercontent.com/ziokim1/QBio305_HHU/main/data/305_linkage_map_data.csv
wget -nc -P ~/Documents/QBIO305_group1_assignment_data https://raw.githubusercontent.com/ziokim1/QBio305_HHU/main/data/305_owb_input_map.csv

echo "Setup Complete"