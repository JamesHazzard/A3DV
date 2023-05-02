#!/bin/bash

copy_config(){

config="./config/jameshome.ini" # choose a config file
cp ${config} scripts/config.ini # copy config file to scripts

}

copy_config
cd scripts
python3 make_sampled_posterior.py

rm -f gmt.* *.cpt