#!/bin/bash

copy_config(){

config="./config/jameshome.ini" # choose a config file
cp ${config} scripts/config.ini # copy config file to scripts

}

copy_config
cd scripts

chmod +x make_divide_grids.sh

./make_divide_grids.sh

rm -f gmt.* junk*