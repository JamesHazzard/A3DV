#!/bin/bash

copy_config(){

config="./config/jameshome.ini" # choose a config file
cp ${config} scripts/config.ini # copy config file to scripts

}

run_make_thermodynamic_outputs(){

    for ((idx=563; idx<1000; idx++)); do
        python3 make_thermodynamic_outputs.py -v distribution -i ${idx}
    done

}

copy_config
cd scripts

chmod +x make_thermodynamic_output_summary.sh
chmod +x make_thermodynamic_output_plots.sh

#run_make_thermodynamic_outputs
#./make_thermodynamic_output_summary.sh
./make_thermodynamic_output_plots.sh

rm -f gmt.* junk*
