gmt gmtset FORMAT_FLOAT_OUT %.9f

################################################################################################################################################
# 1) Make plate model vs. Vs comparison
################################################################################################################################################

make_plate(){

    echo "Making oceanic Vs vs. T file"

    # Remake Vs and Vs uncertainty grids
    maxz=$(gmt grdinfo $Tsgrid -C | awk '{print $5}')
    awk '{print $4, $2, $1}' $Vsmodel | gmt surface -I2.5/10 -T0.25 -G$Vsgrid -R0/${maxt}/12.5/400	#surface at small intervals to smooth out the data (remove noise)
    awk '{print $4, $2, $3}' $Vsmodel | gmt surface -I2.5/10 -T0.25 -G$sdgrid -R0/${maxt}/12.5/400

    # Resample at higher resolution
	rgntot="-R0/${maxt}/12.5/140"
    gmt grdsample $Vsgrid $rgntot -I0.5/0.5 -Gjunk1.grd	#sample the data 'at higher resolution' by interpolating onto a finer mesh
    gmt grdsample $Tsgrid $rgntot -I0.5/0.5 -Gjunk2.grd
    gmt grdsample $sdgrid $rgntot -I0.5/0.5 -Gjunk3.grd

    # Create temporary file with key data
    gmt grd2xyz junk1.grd | awk '{print $0}' > junk1
    gmt grd2xyz junk2.grd | awk '{print $3}' > junk2
    gmt grd2xyz junk3.grd | awk '{print $3}' > junk3
    paste junk1 junk2 junk3 | awk '{print $1, $2, $3, $4, $5}' > junk4

    # Get depth slices
    echo -n > depth_slices
    res=$(echo $zbase_plate $ztop_plate $divs_plate | awk '{print ($1-$2)/$3}')
    echo "Vertical resolution is $res km"
    for i in $(seq $ztop_plate $res $(echo $zbase_plate - $res | bc -l)); do
        zmid=$(echo $i $res | awk '{print $1+($2/2)}')
        maxz=$(echo $i+$res | bc -l)
        echo $zmid $i $maxz >> depth_slices
    done

    # Extract age, depth, Vs, temperature and Vs uncertainty for each depth slice
    echo -n > temps
	echo -n > rest
    depthlines=$(wc -l depth_slices | awk '{print $1}')
    for ((l=1; l<=$depthlines; l++)); do
        for age in $(seq 0 2 $maxt); do
            z=$(awk 'NR=='$l'{printf "%.2f", $1}' depth_slices)
            zmin=$(awk 'NR=='$l'{printf "%.2f", $2}' depth_slices)
            zmax=$(awk 'NR=='$l'{printf "%.2f", $3}' depth_slices)
            awk '{if ($1=='$age' && $2>='$zmin' && $2<='$zmax' && $3!="NaN" && $4!="NaN" && $5!="NaN")\
				{sum+=$3; sum3+=$5;count++}}END{print '$age', '$z', sum/count, sum3/count}' junk4 >> \
				rest
            awk '{if ($1=='$age' && $2=='$z' && $3!="NaN" && $4!="NaN" && $5!="NaN"){print $4}}' junk4 >> \
				temps
        done
    done
	paste rest temps | awk '{print $1, $2, $3, $5, $4}' > ${fold_data_output}/inv_data/plate/plate_${model}${end}.txt

    #  Clean up
    rm -f ages junk* rest temps depth_slices

}

# set misc params
region="global"
poly=${fold_data_input}/exclusion_polygons/mark_polys/oceanic_polygons/${region}.ll
model="ANT-20"
wave="S"
orientation="i"
type="abs"
tomo=${model}_${wave}_${orientation}_${type}
maxt=120
age_inc=2
res=1
Tp_plate=$(awk '$1 ~ /^Tp_plate/' config.ini | awk '{print $3}')
zp_plate=$(awk '$1 ~ /^zp_plate/' config.ini | awk '{print $3}')
rd_plate=$(awk '$1 ~ /^rd_plate/' config.ini | awk '{print $3}')
ztop_plate=$(awk '$1 ~ /^ztop_plate/' config.ini | awk '{print $3}')
zbase_plate=$(awk '$1 ~ /^zbase_plate/' config.ini | awk '{print $3}')
divs_plate=$(awk '$1 ~ /^divs_plate/' config.ini | awk '{print $3}')
end="_avg_noTavg_${ztop_plate}_${zbase_plate}"
z_inc=$(echo $ztop_plate $zbase_plate $divs_plate | awk '{print ($2-$1)/$3}')

echo -n > depth_slices
for ((i=1; i<=$divs_plate; i++)); do
    z=$(echo $ztop_plate $i $z_inc | awk '{print (($2-1)*$3) + $1 + ($3/2)}')
    echo $z >> depth_slices
done

# set up paths to files/dirs
fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_invdata=${fold_base}/inversiondata
fold_scripts=${fold_invdata}/scripts
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_tomo=${fold_base}/${tomo}
f_ocean_age=${fold_data_input}/ocean_age/oc_age_M16_mar18edit.grd
f_ocean_age_sampled=${fold_data_input}/ocean_age/sampled_oc_age_M16_mar18edit.grd
Tsgrid=${fold_data_input}/plate_thermal/Tgrid-${Tp_plate}-${zp_plate}-${rd_plate}.grd
fold_output_velocity=${fold_data_output}/binned_velocity
Vsmodel=${fold_output_velocity}/${region}_${tomo}_mean_velocities_${age_inc}_Ma_bins.vset
Vsgrid=${fold_output_velocity}/${tomo}_vel.grd
sdgrid=${fold_output_velocity}/${tomo}_vel_uncertainty.grd

mkdir -p ${fold_data_output}/inv_data
for output in plate adiabat attenuation viscosity; do
	mkdir -p ${fold_data_output}/inv_data/${output}
done

fold_output_velocity=${fold_data_output}/binned_velocity
mkdir -p ${fold_output_velocity}

make_plate

