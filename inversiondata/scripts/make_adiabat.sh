#!/bin/bash

# code used to prepare adiabatic Vs-T files for use in BANCAL22 inversion

# 1. convert the oceanic age grid into a .dat file, excluding known dodgy regions
# 2. for each location in the tomographic grid, extract location/velocity/age data and save in a .llavz file
# 3. collect all the data within each depth bin and average the tomographic velocity [over age, up to 120 Ma] within each bin

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

extract_velocities () {

# Subroutine to extract velocity data from grids

# Sample the tomography grids at each depth, recording the velocity and depths used
echo -n > ${fold_output_velocity}/${tomo}_depths_used.z
echo "Lon Lat Age Velocity Depth" > ${fold_output_velocity}/${tomo}.llavz
echo "Lon Lat Age Velocity Depth" > ${fold_output_velocity}/${tomo}_fullocean.llavz

for slice in ${fold_tomo}/*km.grd; do

	depth=$(echo ${slice##*/} | sed s/"_km.grd"//g)
    echo "Working on depth slice at "$depth" km..."
    
	filename=$(echo ${slice})
                    
    cat sample_locs.dat | gmt grdtrack -fg -nl -G$filename | grep -v "NaN" | awk '{print $1, $2, $3, $4, '$depth'}' >> ${fold_output_velocity}/${tomo}.llavz
    cat full_locs.dat | gmt grdtrack -fg -nl -G$filename | grep -v "NaN" | awk '{print $1, $2, $3, $4, '$depth'}' >> ${fold_output_velocity}/${tomo}_fullocean.llavz
    echo $depth >> ${fold_output_velocity}/${tomo}_depths_used.z

done

sort -n ${fold_output_velocity}/${tomo}_depths_used.z > a
mv a ${fold_output_velocity}/${tomo}_depths_used.z

}

make_age_bins () {

# Bin the data by age
echo -n > ${Vsmodel}
for age in $(seq 0 $age_inc 200) ; do	#loop over each age bin

	echo "Working on "$age" Ma bin..."
	age_max=$(echo $age $age_inc | awk '{print $1+$2}')
	age_ave=$(echo $age $age_inc | awk '{print $1+($2/2)}')

	# Extract the subset of data in these age bins and also in the basin of interest
	awk 'NR>1 && $3>='$age' && $3<'$age_max' {print $0}' ${fold_output_velocity}/${tomo}.llavz > subset.temp

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

		for depth in $(cat ${fold_output_velocity}/${tomo}_depths_used.z); do
	
			mean=$(awk '$5=='$depth' {print $0}' regional_data.temp | awk '{sum+=$4}END{print sum/NR}')
			sd=$(awk '$5=='$depth' {print $0}' regional_data.temp | awk  '{sum+=$4; sumsq+=$4*$4}END{print sqrt(sumsq/NR - (sum/NR)**2)}')

			echo $mean $depth $sd $age_ave >> ${Vsmodel}	#for each age;depth we store the average (over lon/lat) seismic wave speed and its error
			echo "Determining depth averaged Vs for $depth km in $age age bin of $region region"
	
		done

	fi

	rm subset.temp locs.temp regional_data.temp 2>/dev/null

done

}

make_adiabat(){

    echo "Making adiabat file"

    # Depth range
    min_ad_dep=225
    max_ad_dep=400

    # Extract depth slices
   	ls ${fold_tomo}/ | sed s/"_km.grd"//g | awk '{if ($1<='${max_ad_dep}' && $1>='${min_ad_dep}' && $1%25==0)\
        {printf "%i\n", $1}}' | sort -nk1 > depth_slices

	# Adiabat parameters
	g=9.81
	Cp=1187.
	alpha=3.e-5

    # Extract Vs, depth, temperature and Vs uncertainty for each depth slice
    depthlines=$(wc -l depth_slices | awk '{print $1}')
    echo -n > ${fold_data_output}/inv_data/adiabat/adiabat_${model}.txt
    for ((l=1; l<=$depthlines; l++)); do
	    depth=$(awk 'NR=='$l'{print $1}' depth_slices)
	    Vs=$(awk '{if ($2=='$depth' && $4<='$maxt'){sum=sum+$1; count++}}END{print sum/count}' $Vsmodel)
	    err=$(awk '{if ($2=='$depth' && $4<='$maxt'){sum=sum+$3; count++}}END{print sum/count}' $Vsmodel)
	    z=$(awk 'NR=='$l'{printf "%.1f", $1}' depth_slices)
		T=$(echo $Tp_adiabat $g $Cp $alpha $depth | awk '{print (($1+273)*exp(($2*$4*$5*1000)/$3))-273}')
	    echo $Vs $z $T $err >> ${fold_data_output}/inv_data/adiabat/adiabat_${model}.txt
    done

	#  Clean up
    rm -f depth_slices

}

# set misc params
region="global"
model="ANT-20"
wave="S"
orientation="i"
type="abs"
tomo=${model}_${wave}_${orientation}_${type}
maxt=120
age_inc=2
res=1
Tp_adiabat=$(awk '$1 ~ /^Tp_plate/' config.ini | awk '{print $3}')

# set up paths to files/dirs
fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_invdata=${fold_base}/inversiondata
fold_scripts=${fold_invdata}/scripts
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_tomo=${fold_base}/${tomo}
f_ocean_age=${fold_data_input}/ocean_age/oc_age_M16_mar18edit.grd
f_ocean_age_sampled=${fold_data_input}/ocean_age/sampled_oc_age_M16_mar18edit.grd
fold_output_velocity=${fold_data_output}/binned_velocity
Vsmodel=${fold_output_velocity}/${region}_${tomo}_mean_velocities_${age_inc}_Ma_bins.vset
poly=${fold_data_input}/exclusion_polygons/mark_polys/oceanic_polygons/${region}.ll

mkdir -p ${fold_data_output}/inv_data
for output in plate adiabat attenuation viscosity; do
	mkdir -p ${fold_data_output}/inv_data/${output}
done

mkdir -p ${fold_output_velocity}

get_age_locs
extract_velocities
make_age_bins
make_adiabat