#!/bin/bash

copy_config(){

config="./config/jameshome.ini" # choose a config file
cp ${config} scripts/config.ini # copy config file to scripts

}

run_complex_viscosity_specific_model(){

    for ((idx=0; idx<1000; idx++)); do
        echo "calculating complex viscosity at location ${loc} using stress model ${stress_model} and anelasticity model ${idx}"
        python3 make_relaxation_spectrum.py -v distribution -i ${idx} -l ${loc} -s ${stress_model} -t ${depth_top} -b ${depth_bottom}
    done

}

run_complex_viscosity_all_models(){

    loc="ASE"
    stress_model="B18"
    depth_top="150"
    depth_bottom="175"
    #run_complex_viscosity_specific_model
    
    loc="S20AP_perimeter"
    stress_model="S20"
    depth_top="125"
    depth_bottom="175"
    #run_complex_viscosity_specific_model

    loc="I11AP_perimeter"
    stress_model="I11"
    depth_top="125"
    depth_bottom="175"
    #run_complex_viscosity_specific_model

    loc="Wo15AP_perimeter"
    stress_model="Wo15"
    depth_top="125"
    depth_bottom="175"
    #run_complex_viscosity_specific_model

    loc="S20AP"
    stress_model="S20"
    depth_top="125"
    depth_bottom="175"
    #run_complex_viscosity_specific_model

    loc="I11AP"
    stress_model="I11"
    depth_top="125"
    depth_bottom="175"
    #run_complex_viscosity_specific_model

    loc="Wo15AP"
    stress_model="Wo15"
    depth_top="125"
    depth_bottom="175"
    #run_complex_viscosity_specific_model

    loc="ASE"
    stress_model="B18"
    depth_top="150"
    depth_bottom="250"
    #run_complex_viscosity_specific_model

    loc="S20AP_perimeter"
    stress_model="S20"
    depth_top="125"
    depth_bottom="250"
    #run_complex_viscosity_specific_model

    loc="I11AP_perimeter"
    stress_model="I11"
    depth_top="125"
    depth_bottom="250"
    #run_complex_viscosity_specific_model

    loc="Wo15AP_perimeter"
    stress_model="Wo15"
    depth_top="125"
    depth_bottom="250"
    #run_complex_viscosity_specific_model

    loc="S20AP"
    stress_model="S20"
    depth_top="125"
    depth_bottom="250"
    #run_complex_viscosity_specific_model

    loc="I11AP"
    stress_model="I11"
    depth_top="125"
    depth_bottom="250"
    #run_complex_viscosity_specific_model

    loc="Wo15AP"
    stress_model="Wo15"
    depth_top="125"
    depth_bottom="250"
    #run_complex_viscosity_specific_model

}

run_complex_viscosity_summary_specific_model(){

    python3 make_viscosity_summary.py -l ${loc} -s ${stress_model} -t ${depth_top} -b ${depth_bottom}

}

run_complex_viscosity_summary_all_models(){

    loc="ASE"
    stress_model="B18"
    depth_top="150"
    depth_bottom="175"
    run_complex_viscosity_summary_specific_model
    
    loc="S20AP_perimeter"
    stress_model="S20"
    depth_top="125"
    depth_bottom="175"
    run_complex_viscosity_summary_specific_model

    loc="I11AP_perimeter"
    stress_model="I11"
    depth_top="125"
    depth_bottom="175"
    run_complex_viscosity_summary_specific_model

    loc="Wo15AP_perimeter"
    stress_model="Wo15"
    depth_top="125"
    depth_bottom="175"
    run_complex_viscosity_summary_specific_model

    loc="S20AP"
    stress_model="S20"
    depth_top="125"
    depth_bottom="175"
    run_complex_viscosity_summary_specific_model

    loc="I11AP"
    stress_model="I11"
    depth_top="125"
    depth_bottom="175"
    run_complex_viscosity_summary_specific_model

    loc="Wo15AP"
    stress_model="Wo15"
    depth_top="125"
    depth_bottom="175"
    run_complex_viscosity_summary_specific_model

    loc="ASE"
    stress_model="B18"
    depth_top="150"
    depth_bottom="250"
    run_complex_viscosity_summary_specific_model

    loc="S20AP_perimeter"
    stress_model="S20"
    depth_top="125"
    depth_bottom="250"
    run_complex_viscosity_summary_specific_model

    loc="I11AP_perimeter"
    stress_model="I11"
    depth_top="125"
    depth_bottom="250"
    run_complex_viscosity_summary_specific_model

    loc="Wo15AP_perimeter"
    stress_model="Wo15"
    depth_top="125"
    depth_bottom="250"
    run_complex_viscosity_summary_specific_model

    loc="S20AP"
    stress_model="S20"
    depth_top="125"
    depth_bottom="250"
    run_complex_viscosity_summary_specific_model

    loc="I11AP"
    stress_model="I11"
    depth_top="125"
    depth_bottom="250"
    run_complex_viscosity_summary_specific_model

    loc="Wo15AP"
    stress_model="Wo15"
    depth_top="125"
    depth_bottom="250"
    run_complex_viscosity_summary_specific_model

}

copy_config
cd scripts

#chmod +x make_sample_loc_grids.sh
#chmod +x make_sample_loc_plots.sh
#chmod +x make_viscosity_summary_plots.sh

#./make_sample_loc_grids.sh
#./make_sample_loc_plots.sh
#run_complex_viscosity_all_models
#run_complex_viscosity_summary_all_models
./make_viscosity_summary_plots.sh

#loc="S21AP"
#stress_model="S21"
#depth_top="125"
#depth_bottom="175"
#run_complex_viscosity_specific_model
#run_complex_viscosity_summary_specific_model

rm -f gmt.* junk* *.cpt
