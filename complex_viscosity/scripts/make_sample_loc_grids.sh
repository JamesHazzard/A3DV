#!/bin/bash

sample_temperature_grids_by_perimeter(){

    # use gmt select to extract regions within lon/lat
    f_sample_locs=${fold_data_input}/${input_locs}/sample_perimeter_locs.txt
    no_locs=$(wc -l ${f_sample_locs} | awk '{print $1}')

    #for ((i=1; i<=1; i++)); do
    for ((i=1; i<=${no_locs}; i++)); do

        j=$(echo $i | awk '{print $1+1}')
        loc=$(awk '{if (NR=='${j}') print $1}' ${f_sample_locs})
        loc_full=$(awk '{if (NR=='${j}') print $2}' ${f_sample_locs})

        fold_data_output_loc=${fold_data_output}/location_${loc_full}
        mkdir -p ${fold_data_output_loc}

        n=$(wc -l ${fold_data_input}/${input_locs}/${loc}_perimeter.txt | awk '{print $1-1}')
        perimeter="junk_polygon.txt"
        tail -n ${n} ${fold_data_input}/${input_locs}/${loc}_perimeter.txt > ${perimeter}

        #for depth in $(seq 100 25 100); do
        for depth in $(seq 100 25 250); do

            check_grid_name=depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_0
            f_check_grid=${fold_data_input_temperature}/distribution/individual/${check_grid_name}.grd

            if [[ -e $f_check_grid ]]; then

                for ((k=0; k<1000; k++)); do

                    echo "sampling temperature grid ${k} at location ${loc} for depth slice ${depth}.."
                
                    idx=$k
                    grid_name=depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_${idx}
                    f_grid=${fold_data_input_temperature}/distribution/individual/${grid_name}.grd
                    out_name=${fold_data_output_loc}/${grid_name}_${loc_full}.xyz
                    if [[ ! -e ${out_name} ]]; then
                        gmt grd2xyz -R285/307/-77/-61 ${f_grid} | gmt gmtselect -F${perimeter} -fg > ${out_name}
                        #gmt gmtselect ${f_grid} -F${perimeter} -bi -fg > ${out_name}
                    fi 

                done

            fi
                
        done

    done

}

sample_viscosity_grids(){

    fold_data_output=${fold_data_output_viscosity}/${date}/Tp_${Tp}_sol50_${sol50}
    f_sample_locs=${fold_data_input}/${input_locs}/sample_locs.txt
    no_locs=$(wc -l ${f_sample_locs} | awk '{print $1}')
    rad=50k
    
    # loop over each location we wish to sample, extracting viscosity as a function of space and anelasticity model
    for ((i=1; i<=${no_locs}; i++)); do
    #for ((i=1; i<=1; i++)); do

        j=$(echo $i | awk '{print $1+1}')
        loc=$(awk '{if (NR=='${j}') print $1}' ${f_sample_locs})
        lon=$(awk '{if (NR=='${j}') print $2}' ${f_sample_locs})
        lat=$(awk '{if (NR=='${j}') print $3}' ${f_sample_locs})

        fold_data_output_loc=${fold_data_output}/location_${loc}
        mkdir -p ${fold_data_output_loc}

        for depth in $(seq 100 25 400); do

            check_grid_name=depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_0
            f_check_grid=${fold_data_input_viscosity}/distribution/individual/${check_grid_name}.grd

            if [[ -e $f_check_grid ]]; then

                echo "sampling viscosity grids at location ${loc} for depth slice ${depth}.."

                for ((k=0; k<1000; k++)); do
                
                    idx=$k
                    grid_name=depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_${idx}
                    f_grid=${fold_data_input_viscosity}/distribution/individual/${grid_name}.grd
                    gmt grdcut -S${lon}/${lat}/${rad} ${f_grid} -G${fold_data_output_loc}/${grid_name}_${loc}.grd
                    #gmt grd2xyz ${fold_data_output_loc}/${grid_name}_${loc}.grd > ${fold_data_output_loc}/${grid_name}_${loc}.xyz

                done

            fi
            
        done

    done

    echo "done"

}

sample_temperature_grids(){

    fold_data_output=${fold_data_output_temperature}/${date}/Tp_${Tp}_sol50_${sol50}
    f_sample_locs=${fold_data_input}/${input_locs}/sample_locs.txt
    no_locs=$(wc -l ${f_sample_locs} | awk '{print $1}')
    rad=50k
    
    # loop over each location we wish to sample, extracting temperature as a function of space and anelasticity model
    for ((i=4; i<=4; i++)); do
    #for ((i=1; i<=${no_locs}; i++)); do

        j=$(echo $i | awk '{print $1+1}')
        loc=$(awk '{if (NR=='${j}') print $1}' ${f_sample_locs})
        lon=$(awk '{if (NR=='${j}') print $2}' ${f_sample_locs})
        lat=$(awk '{if (NR=='${j}') print $3}' ${f_sample_locs})

        fold_data_output_loc=${fold_data_output}/location_${loc}
        mkdir -p ${fold_data_output_loc}

        #for depth in 275 300; do
        for depth in $(seq 100 25 250); do

            check_grid_name=depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_0
            f_check_grid=${fold_data_input_temperature}/distribution/individual/${check_grid_name}.grd

            if [[ -e $f_check_grid ]]; then

                echo "sampling temperature grids at location ${loc} for depth slice ${depth}.."

                for ((k=0; k<1000; k++)); do
                
                    idx=$k
                    grid_name=depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_${idx}
                    f_grid=${fold_data_input_temperature}/distribution/individual/${grid_name}.grd
                    out_name=${fold_data_output_loc}/${grid_name}_${loc}.grd
                    if [[ ! -e ${out_name} ]]; then
                        gmt grdcut -S${lon}/${lat}/${rad} ${f_grid} -G${fold_data_output_loc}/${grid_name}_${loc}.grd
                        gmt grd2xyz ${fold_data_output_loc}/${grid_name}_${loc}.grd > ${fold_data_output_loc}/${grid_name}_${loc}.xyz
                    fi 

                done

            fi
            
        done

    done

    echo "done"

}

fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
input_locs=$(awk '$1 ~ /^input_locations/' config.ini | awk '{print $3}')
input_thermodynamic=$(awk '$1 ~ /^input_thermodynamic/' config.ini | awk '{print $3}')
input_viscosity=$(awk '$1 ~ /^input_viscosity/' config.ini | awk '{print $3}')
input_temperature=$(awk '$1 ~ /^input_temperature/' config.ini | awk '{print $3}')
fold_data_input_thermodynamic=${fold_data_input}/${input_thermodynamic}
# since this is a data preparation stage, data output will be routed towards an input dir
fold_data_output_viscosity=${fold_data_input}/${input_viscosity}
fold_data_output_temperature=${fold_data_input}/${input_temperature}
Tp=$(awk '$1 ~ /^potential_temperature/' config.ini | awk '{print $3}')
sol50=$(awk '$1 ~ /^solidus_50km/' config.ini | awk '{print $3}')
date=$(awk '$1 ~ /^date/' config.ini | awk '{print $3}')
fold_data_input_viscosity=${fold_data_input_thermodynamic}/${date}/Tp_${Tp}_sol50_${sol50}/viscosity_grids
fold_data_input_temperature=${fold_data_input_thermodynamic}/${date}/Tp_${Tp}_sol50_${sol50}/temperature_grids

fold_data_output=${fold_data_output_viscosity}
mkdir -p ${fold_data_output}
fold_data_output=${fold_data_output}/${date}
mkdir -p ${fold_data_output}
fold_data_output=${fold_data_output}/Tp_${Tp}_sol50_${sol50}
mkdir -p ${fold_data_output}

fold_data_output=${fold_data_output_temperature}
mkdir -p ${fold_data_output}
fold_data_output=${fold_data_output}/${date}
mkdir -p ${fold_data_output}
fold_data_output=${fold_data_output}/Tp_${Tp}_sol50_${sol50}
mkdir -p ${fold_data_output}

sample_temperature_grids
#sample_temperature_grids_by_perimeter
#sample_viscosity_grids