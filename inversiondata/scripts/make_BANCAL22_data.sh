#!/bin/bash 

remove_existing_data(){

    # deletes current inversion dir data set for fresh transfer (important if adding/subtracting data types e.g. removing viscosity constraint for an inv run)

    cd ${BANCALfold}

    for data_set in plate adiabat attenuation viscosity; do
        cd ${data_set}
        rm -f *
        cd ${BANCALfold}
    done

}

make_data(){

    # formats data ready for inversion

    model="ANT-20"
    Qmodel="QRFSI12"
    end="avg_noTavg_50_125"

    awk '{print $3, $5, $4, $2}' ${fold_data_output}/plate/plate_${model}_${end}.txt > ${fold_data_output}/BANCAL22_data/plate.VseTz
    awk '{print $1, $4, $3, $2}' ${fold_data_output}/adiabat/adiabat_${model}.txt > ${fold_data_output}/BANCAL22_data/adiabat.VseTz
    awk '{print $1, $5, $2, $3}' ${fold_data_output}/attenuation/attenuation_${Qmodel}.txt > ${fold_data_output}/BANCAL22_data/attenuation.QeVsz
    awk '{print $1, $5, $2, $3}' ${fold_data_output}/viscosity/viscosity_${model}.txt > ${fold_data_output}/BANCAL22_data/viscosity.neVsz
    echo ${Tp} > ${fold_data_output}/BANCAL22_data/Tp_${model}.txt
    echo ${sol50} > ${fold_data_output}/BANCAL22_data/sol50_${model}.txt

}

transfer_data(){

    # transfers data set to inversion dir for running

    model="ANT-20"
    cp ${fold_data_output}/BANCAL22_data/plate.VseTz ${BANCALfold}/plate/plate.VseTz
    cp ${fold_data_output}/BANCAL22_data/adiabat.VseTz ${BANCALfold}/adiabat/adiabat.VseTz
    cp ${fold_data_output}/BANCAL22_data/attenuation.QeVsz ${BANCALfold}/attenuation/attenuation.QeVsz
    cp ${fold_data_output}/BANCAL22_data/viscosity.neVsz ${BANCALfold}/viscosity/viscosity.neVsz
    cp ${fold_data_output}/BANCAL22_data/Tp_${model}.txt ${BANCALfold}/potential_temperature/potential_temperature.T
    cp ${fold_data_output}/BANCAL22_data/sol50_${model}.txt ${BANCALfold}/potential_temperature/solidus_50km_temperature.T
}

# set misc params
Tp=$(awk '$1 ~ /^Tp_plate/' config.ini | awk '{print $3}')
sol50=$(awk '$1 ~ /^solidus_50km/' config.ini | awk '{print $3}')

# set up paths to files/dirs
fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_invdata=${fold_base}/inversiondata
fold_scripts=${fold_invdata}/scripts
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')/inv_data
BANCALfold=$(awk '$1 ~ /^BANCAL22/' config.ini | awk '{print $3}')

mkdir -p ${fold_data_output}/BANCAL22_data

make_data


