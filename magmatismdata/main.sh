#!/bin/bash

copy_config(){

config="./config/jameshome.ini" # choose a config file
cp ${config} scripts/config.ini # copy config file to scripts

}

copy_config
cd scripts

chmod +x make_LAB_magmatism_comparison.sh
chmod +x make_plot_LAB_magmatism.sh

echo "extracting magmatism data.."
#python3 extract_data.py
echo "done"
echo "extracting relationship between LAB depth and age since magmatism.."
#./make_LAB_magmatism_comparison.sh
echo "done"
echo "calculating correlation coefficient by simulating range of data sets.. (this may take ~10 mins if running for first time)"
#python3 calculate_spearman.py
echo "done"
echo "plotting results.."
./make_plot_LAB_magmatism.sh
echo "done"

rm -f gmt.* junk*