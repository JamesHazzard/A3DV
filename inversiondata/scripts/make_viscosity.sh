#!/bin/bash

make_viscosity(){

    echo "Making viscosity file"

	#Adiabat parameters
	g=9.81
	Cp=1187.
	alpha=3.e-5

    min_visc_dep=220
    max_visc_dep=400

    maxt=120

    # Extract depth slices
   	ls ${fold_tomo}/ | sed s/"_km.grd"//g | awk '{if ($1<='${max_visc_dep}' && $1>='${min_visc_dep}' && $1%25==0)\
        {printf "%i\n", $1}}' | sort -nk1 > depth_slices

    # Extract viscosity, Vs, depth, temperature, viscosity uncertainty and Vs uncertainty for each depth slice
    shallow_eta=$(echo 1.e20 | awk '{print log($1)/log(10)}')
    depthlines=$(wc -l depth_slices | awk '{print $1}')
    echo -n > ${fold_data_output}/inv_data/viscosity/viscosity_${model}.txt
    for ((l=1; l<=$depthlines; l++)); do
	    depth=$(awk 'NR=='$l'{print $1}' depth_slices)
	    Vs=$(awk '{if ($2=='$depth' && $4<='$maxt'){sum=sum+$1; count++}}END{print sum/count}' $Vsmodel)
	    err=$(awk '{if ($2=='$depth' && $4<='$maxt'){sum=sum+$3; count++}}END{print sum/count}' $Vsmodel)
	    z=$(awk 'NR=='$l'{printf "%.1f", $1}' depth_slices)
		T=$(echo $Tp_viscosity $g $Cp $alpha $z | awk '{print (($1+273)*exp(($2*$4*$5*1000)/$3))-273}')
	    echo $shallow_eta $Vs $z $T "1" $err >> ${fold_data_output}/inv_data/viscosity/viscosity_${model}.txt
    done

    #  Clean up
    rm depth_slices

}

# Set misc params
region="global"
age_inc=2
model="ANT-20"
wave="S"
orientation="i"
type="abs"
tomo="${model}_${wave}_${orientation}_${type}"
Tp_viscosity=$(awk '$1 ~ /^Tp_plate/' config.ini | awk '{print $3}')

# set up paths to files/dirs
fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_invdata=${fold_base}/inversiondata
fold_scripts=${fold_invdata}/scripts
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_tomo=${fold_base}/${tomo}
fold_output_velocity=${fold_data_output}/binned_velocity
Vsmodel=${fold_output_velocity}/${region}_${tomo}_mean_velocities_${age_inc}_Ma_bins.vset

make_viscosity

