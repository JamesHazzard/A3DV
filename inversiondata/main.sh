#!/bin/bash

copy_config(){

config="./config/jameshome.ini" # choose a config file
cp ${config} scripts/config.ini # copy config file to scripts

}

copy_config
cd scripts

chmod +x make_adiabat.sh
chmod +x make_plate.sh
chmod +x make_attenuation.sh
chmod +x make_viscosity.sh
chmod +x make_BANCAL22_data.sh

./make_adiabat.sh
./make_plate.sh
./make_attenuation.sh
./make_viscosity.sh
./make_BANCAL22_data.sh

rm -f *_locs.dat subset.temp gmt.* *.temp