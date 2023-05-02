#!/bin/bash

copy_config(){

config="./config/jameshome.ini" # choose a config file
cp ${config} scripts/config.ini # copy config file to scripts

}

copy_config
cd scripts  # change working directory to location of scripts

fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}') # create log file
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')

chmod +x make_plot_pappa_comparison.sh

./make_plot_pappa_comparison.sh

rm -f gmt. *.cpt *junk*
