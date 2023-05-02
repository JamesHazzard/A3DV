#!/bin/bash

copy_config(){

config="./config/jameshome.ini" # choose a config file
cp ${config} scripts/config.ini # copy config file to scripts

}

copy_config
cd scripts

#python3 make_parameter_summary.py
#chmod +x make_trade_off_density_plot.sh
#./make_trade_off_density_plot.sh
#chmod +x make_parameter_histograms.sh
#./make_parameter_histograms.sh
#chmod +x make_credible_interval_data_fits.sh
#./make_credible_interval_data_fits.sh
#python3 make_credible_interval_envelopes.py
#chmod +x make_credible_interval_data_fits_plot.sh
./make_credible_interval_data_fits_plot.sh

wait
rm -f junk* gmt.* *.cpt *depth_slices*