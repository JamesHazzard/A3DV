#!/bin/bash

get_age_locs () {

# Subroutine to find locations in age grid at which Vs will be sampled
# Downsize the grid and mask points inside the dodgy polygons
gmt grdsample ${f_ocean_age} -I$res -G${f_ocean_age_sampled}
# Convert to sample points and clip the regions of known dodgy oceanic crust
echo -n > polys.temp
# For cut_out in /space1/mjh217/residual_depth/shiptrack/anomalous_crust/mark_polys/*.mask ; do
for cut_out in ${fold_data_input}/exclusion_polygons/mark_polys/*.mask ; do
	echo ">" >> polys.temp
	cat $cut_out >> polys.temp
done
# NB. Double set of nodes at 180 and -180 that don't get clipped by exclusion polys. Just throw out now (enough other good stuff).
gmt grd2xyz ${f_ocean_age_sampled} -fg | grep -v "NaN" | awk '$1>-180 && $1<180 {print $1, $2, $3}' | gmt gmtselect -fg -Fpolys.temp -If > sample_locs.dat 
gmt grd2xyz ${f_ocean_age_sampled} -fg | grep -v "NaN" | awk '$1>-180 && $1<180 {print $1, $2, $3}' > full_locs.dat 
rm polys.temp

}

extract_Q () {

# Subroutine to extract attenuation data from grids

# Sample the tomography grids at each dewpth, recording the attenuation and depths used
echo -n > ${fold_output_attenuation}/${Qmod}_depths_used.z
echo "Lon Lat Age Attenuation Depth" > ${fold_output_attenuation}/${Qmod}.llaQz
echo "Lon Lat Age Attenuation Depth" > ${fold_output_attenuation}/${Qmod}_fullocean.llaQz

for slice in ${fold_Qmod}/*km.grd; do

	depth=$(echo ${slice##*/} | sed s/"_km.grd"//g)
    echo "Working on depth slice at "$depth" km..."

	filename=$(echo ${slice})
                    
    cat sample_locs.dat | gmt grdtrack -fg -nl -G$filename | grep -v "NaN" | awk '{print $1, $2, $3, $4, '$depth'}' >> ${fold_output_attenuation}/${Qmod}.llaQz
    cat full_locs.dat | gmt grdtrack -fg -nl -G$filename | grep -v "NaN" | awk '{print $1, $2, $3, $4, '$depth'}' >> ${fold_output_attenuation}/${Qmod}_fullocean.llaQz
    echo $depth >> ${fold_output_attenuation}/${Qmod}_depths_used.z

done

sort -n ${fold_output_attenuation}/${Qmod}_depths_used.z > a
mv a ${fold_output_attenuation}/${Qmod}_depths_used.z

}

make_age_bins () {

# Bin the data by age
echo -n > ${Qmodel}
for age in $(seq 0 $age_inc 200) ; do	#loop over each age bin

	echo "Working on "$age" Ma bin..."
	age_max=$(echo $age $age_inc | awk '{print $1+$2}')
	age_ave=$(echo $age $age_inc | awk '{print $1+($2/2)}')

	# Extract the subset of data in these age bins and also in the basin of interest
	awk 'NR>1 && $3>='$age' && $3<'$age_max' {print $0}' ${fold_output_attenuation}/${Qmod}.llaQz > subset.temp

	if [ $region == "global" ] ; then
		awk '{print $1, $2}' subset.temp | sort -nk1 -nk2 | uniq | gmt gmtselect \
      -F${fold_data_input}/exclusion_polygons/mark_polys/oceanic_polygons/atlantic.ll -fg > locs.temp
		awk '{print $1, $2}' subset.temp | sort -nk1 -nk2 | uniq | gmt gmtselect \
		  -F${fold_data_input}/exclusion_polygons/mark_polys/oceanic_polygons/pacific.ll -fg >> locs.temp
		awk '{print $1, $2}' subset.temp | sort -nk1 -nk2 | uniq | gmt gmtselect \
		  -F${fold_data_input}/exclusion_polygons/mark_polys/oceanic_polygons/indian.ll -fg >> locs.temp
	else
		awk '{print $1, $2}' subset.temp | sort -nk1 -nk2 | uniq | gmt gmtselect -F$poly -fg > locs.temp
	fi

	# Loop over the unique locations

	n=$(wc -l locs.temp | awk '{print $1}')

	if [ $n -eq 0 ] ; then

		echo "No points in age range "$age" to "$age_max" Ma"

	else

		echo -n > regional_data.temp
		for ((i=1 ; i<=$n; i++)); do	#loop over each location inside each age bin 

			lon=$(awk 'NR=='$i' {print $1}' locs.temp)
			lat=$(awk 'NR=='$i' {print $2}' locs.temp)

			awk '$1=='$lon' && $2=='$lat' {print $0}' subset.temp >> regional_data.temp	#extract the location of each data point within each age;depth bin 
			echo "Extracting data point $i of $n in $age Ma age bin of $region region"

		done

		# Plot the mean and standard deviation at each depth slice

		for depth in $(cat ${fold_output_attenuation}/${Qmod}_depths_used.z); do
	
			mean=$(awk '$5=='$depth' {print $0}' regional_data.temp | awk '{sum+=$4}END{print sum/NR}')
			sd=$(awk '$5=='$depth' {print $0}' regional_data.temp | awk  '{sum+=$4; sumsq+=$4*$4}END{print sqrt(sumsq/NR - (sum/NR)**2)}')

			echo $mean $depth $sd $age_ave >> ${Qmodel}	#for each age;depth we store the average (over lon/lat) seismic wave speed and its error
			echo "Determining depth averaged Q for $depth km in $age age bin of $region region"
	
		done

	fi

	rm subset.temp locs.temp regional_data.temp 2>/dev/null

done

}

make_attenuation(){

    echo "Making attenuation file"

    # Remake Q^-1 and Q^-1 uncertainty grids
    awk '{print $4, $2, $1}' $Qmodel | gmt surface -I2.5/10 -T0.25 -G$Qgrid -R0/${maxt}/12.5/400
    awk '{print $4, $2, $3}' $Qmodel | gmt surface -I2.5/10 -T0.25 -G$Qsdgrid -R0/${maxt}/12.5/400

    # Resample at higher resolution
	rgntot="-R0/${maxt}/12.5/400"
    gmt grdsample $Qgrid $rgntot -I0.5/0.5 -Gjunk1.grd
    gmt grdsample $Qsdgrid $rgntot -I0.5/0.5 -Gjunk2.grd

    # Create temporary file with key data
    gmt grd2xyz junk1.grd | awk '{print $0}' > junk1
    gmt grd2xyz junk2.grd | awk '{print $3}' > junk2
    paste junk1 junk2 | awk '{print $1, $2, $3, $4}' > junk3

    # Extract depth slices
    min_Q_dep=150
    max_Q_dep=400
	interval_Q_dep=25
	n_intervals=$(echo $min_Q_dep $max_Q_dep $interval_Q_dep | awk '{print ($2-$1)/$3}')
	echo -n > depth_slices
	for ((j=0; j<=$n_intervals; j++)); do
		z=$(echo $min_Q_dep $interval_Q_dep $j | awk '{print $1+($2*$3)}')
		echo $z >> depth_slices
	done

    # Extract age, depth, attenuation, temperature and attenuation uncertainty for each depth slice
    echo -n > attenuation.temp
    agemin=100
    depthlines=$(wc -l depth_slices | awk '{print $1}')
    for ((l=1; l<=$depthlines; l++)); do
        z=$(awk 'NR=='$l'{printf "%.2f", $1}' depth_slices)
        awk '{if ($1>='$agemin' && $1<='$maxt' && $2=='$z' && $3!="NaN" && $4!="NaN")\
            {sum+=$3; sum3+=$4;count++}}END{print '$z', sum/count, sum3/count}' junk3 >> attenuation.temp
    done


	#Adiabat parameters
	g=9.81
	Cp=1187.
	alpha=3.e-5

    # Extract attenuation, Vs, depth, temperature, attenuation uncertainty and Vs uncertainty for each depth slice
    echo -n > ${fold_data_output}/inv_data/attenuation/attenuation_${Qmod}.txt
    lines=$(wc -l depth_slices | awk '{print $1}')
    for ((l=1; l<=$lines; l++)); do
    	z=$(awk 'NR=='$l'{print $1}' depth_slices)
		T=$(echo $Tp_attenuation $g $Cp $alpha $z | awk '{print (($1+273)*exp(($2*$4*$5*1000)/$3))-273}')
   	 	Vs=$(awk '{if ($4>=100 && $4<='$maxt' && $2=='$z'){sum=sum+$1; count++}}END{print sum/count}' $Vsmodel)
   	 	err=$(awk '{if ($4>=100 && $2=='$z'){sum=sum+$3; count++}}END{print sum/count}' $Vsmodel)
    	awk '{if ($1=='$z'){print $2, '$Vs', $1, '$T', $3, '$err'}}' attenuation.temp >> \
			${fold_data_output}/inv_data/attenuation/attenuation_${Qmod}.txt
    done

    #  Clean up
    rm -f depth_slices *.temp junk*

}

# set misc params
region="global"
maxt=120
age_inc=2
res=1
Tp_attenuation=$(awk '$1 ~ /^Tp_plate/' config.ini | awk '{print $3}')
# tomography model params
model="ANT-20"
wave="S"
orientation="i"
type="abs"
tomo=${model}_${wave}_${orientation}_${type}
# attenuation model params
Qmodname="QRFSI12"
Qtype="abs"
Qmod=${Qmodname}
Qbottomz=800
Qregion="global"

# set up paths to files/dirs
fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_invdata=${fold_base}/inversiondata
fold_scripts=${fold_invdata}/scripts
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_tomo=${fold_base}/${tomo}
fold_Qmod=${fold_base}/${Qmod}
f_ocean_age=${fold_data_input}/ocean_age/oc_age_M16_mar18edit.grd
f_ocean_age_sampled=${fold_data_input}/ocean_age/sampled_oc_age_M16_mar18edit.grd
fold_output_velocity=${fold_data_output}/binned_velocity
fold_output_attenuation=${fold_data_output}/binned_attenuation
Vsmodel=${fold_output_velocity}/${region}_${tomo}_mean_velocities_${age_inc}_Ma_bins.vset
Qmodel=${fold_output_attenuation}/${Qregion}_${Qmod}_mean_Q_${age_inc}_Ma_bins.Qet
Qgrid=${fold_output_attenuation}/${Qmod}_attenuation.grd
Qsdgrid=${fold_output_attenuation}/${Qmod}_attenuation_uncertainty.grd
poly=${fold_data_input}/exclusion_polygons/mark_polys/oceanic_polygons/${region}.ll

mkdir -p ${fold_data_output}/inv_data
for output in plate adiabat attenuation viscosity; do
	mkdir -p ${fold_data_output}/inv_data/${output}
done

mkdir -p ${fold_output_velocity}
mkdir -p ${fold_output_attenuation}

get_age_locs
extract_Q
make_age_bins
make_attenuation