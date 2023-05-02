#!/bin/bash

process_magmatism_data(){

  resolution_blockmean=$(awk '{print $1}' ${fold_input_LAB}/${model}/resolution.xx)
  resolution_select=$(awk '{print $2}' ${fold_input_LAB}/${model}/resolution.xx)
  # First extract the location and age information from the two geochemical data sets
  awk -F, '{if(NR>1 && $5<65) print ($3 + $4)/2, ($1 + $2)/2, $5, $6}' ${fold_data_output}/ANTARCTICA_rift_filtered.csv > ${fold_data_output}/Antarctica_magmatism_GEOROC.llas
  awk -F, '{if(NR>1 && $5<65) print $2, $1, $3, $4}' ${fold_data_output}/BALL_filtered.csv > ${fold_data_output}/Antarctica_magmatism_BALL.llas
  awk '{print $1, $2, $3, $4}' ${fold_data_output}/Antarctica_magmatism_GEOROC.llas > ${fold_data_output}/Antarctica_magmatism_aggregated.llas
  awk '{print $1, $2, $3, $4}' ${fold_data_output}/Antarctica_magmatism_BALL.llas >> ${fold_data_output}/Antarctica_magmatism_aggregated.llas
  length=$(awk '{print $0}' ${fold_data_output}/Antarctica_magmatism_aggregated.llas | wc -l)
  echo -n > ints
  for ((i=1; i<=$length; i++)); do
    echo $i >> ints
  done
  paste ${fold_data_output}/Antarctica_magmatism_aggregated.llas ints > ${fold_data_output}/Antarctica_magmatism_aggregated.llasi
  # Now calculate the minimum time since last eruption in each spatial bin (block) - outputs mean location
  awk '{print $1, $2, $3}' ${fold_data_output}/Antarctica_magmatism_aggregated.llas |\
   gmt blockmean -I${resolution_blockmean}k -R-180/180/-90/-60 -E > ${fold_data_output}/Antarctica_magmatism_agg_blockmean.xyzslh
  # Record location of minimum magmatism value in each block and assctd. record ID
  awk '{print $1, $2, $3, $5}' ${fold_data_output}/Antarctica_magmatism_aggregated.llasi |\
   gmt blockmedian -T0.000000000000001 -I${resolution_blockmean}k -R-180/180/-90/-60 -E -Es | awk '{print $1, $2, $4}' > ${fold_data_output}/blockmedian_locs_magmatism.xy
  # Record location of center of each block
  awk '{print $1, $2, $3}' ${fold_data_output}/Antarctica_magmatism_aggregated.llas |\
   gmt blockmean -I${resolution_blockmean}k -R-180/180/-90/-60 -C | awk '{print $1, $2}' > ${fold_data_output}/blockmean_locs_magmatism.xy
  # Create ${fold_data_output}/pointfile using mean locations of magmatism data
  awk '{print $1, $2}' ${fold_data_output}/Antarctica_magmatism_agg_blockmean.xyzslh > ${fold_data_output}/pointfile
  # Select LAB data within 100km radius of mean magmatic locations, and calculate average/std, recording central block location
  gmt select ${LAB_mean}.xyz -C${fold_data_output}/pointfile+d${resolution_select}k -fg | gmt blockmean -I${resolution_blockmean}k -R-180/180/-90/-60 -E -C > ${fold_data_output}/LAB_${model}_blockmean.xyzslh
  gmt select ${LAB_std}.xyz -C${fold_data_output}/pointfile+d${resolution_select}k -fg | gmt blockmean -I${resolution_blockmean}k -R-180/180/-90/-60 -E -C > ${fold_data_output}/LAB_${model}_std_blockmean.xyzslh
  # Match the magmatic data (avg location, min time) back up with the LAB data (avg/std) using the central block location as the connection
  no_locs=$(awk '{print $0}' ${fold_data_output}/blockmean_locs_magmatism.xy | wc -l)
  echo -n > ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean.xyaszs
  echo -n > ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean_range.xyaszs

  for ((loc=1; loc<=$no_locs; loc++)); do

    # Find the longitude and latitude corresponding to the [loc]th data point (location of center of block)
    lon=$(awk '{if(NR=='$loc') print $1}' ${fold_data_output}/blockmean_locs_magmatism.xy)
    lat=$(awk '{if(NR=='$loc') print $2}' ${fold_data_output}/blockmean_locs_magmatism.xy)
    # For the [loc]th data point, extract the mean magmatic location, and the minimum time since last eruption
    awk '{if(NR=='$loc') print $1, $2, $5}' ${fold_data_output}/Antarctica_magmatism_agg_blockmean.xyzslh > junk1
    # Also extract the age uncertainty assctd. with the minimum time since last eruption
    id=$(awk '{if(NR=='$loc') print $3}' ${fold_data_output}/blockmedian_locs_magmatism.xy)
    awk '{if(NR=='$id') print $4}' ${fold_data_output}/Antarctica_magmatism_aggregated.llas > junk2
    # Also extract the mean and std LAB1200 depth
    awk '{if($1=='$lon' && $2=='$lat') print '$loc', $3, $4}' ${fold_data_output}/LAB_${model}_blockmean.xyzslh
    awk '{if($1=='$lon' && $2=='$lat') print $3, $4}' ${fold_data_output}/LAB_${model}_blockmean.xyzslh > junk3
    # Also extract the mean LAB1200 depth uncertainty
    awk '{if($1=='$lon' && $2=='$lat') print $3}' ${fold_data_output}/LAB_${model}_std_blockmean.xyzslh > junk4
    # Add the blockmean std and anelasticity-based std in quadrature
    paste junk3 junk4 | awk '{print $1, ($2^2 + $3^2)^(0.5)}' > junk5
    # Output mean location, minimum age since last eruption, age uncertainty, mean LAB1200 depth, LAB1200 depth uncertainty (spatial + physics)
    paste junk1 junk2 junk5 >> ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean.xyaszs
    # Alternatively use range of location based LAB1200 depth rather than std 
    rm junk3 junk5
    awk '{if($1=='$lon' && $2=='$lat') print $3, $5, $6}' ${fold_data_output}/LAB_${model}_blockmean.xyzslh > junk3
    paste junk3 junk4 | awk '{print $1, ((0.5*($3 - $2))^2 + $4^2)^(0.5)}' > junk5 
    paste junk1 junk2 junk5 >> ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean_range.xyaszs

  done

  rm -f junk* ${fold_data_output}/pointfile ${fold_data_output}/*.llas ${fold_data_output}/LAB*.xyzslh ${fold_data_output}/*locs* ${fold_data_output}/ints ${fold_data_output}/*.xyzslh ${fold_data_output}/*.llasi

}

plot_LAB_magmatism_map(){

  ps="./plots/LAB1200_${model}_magmatism_map.ps"
  jpg="./plots/LAB1200_${model}_magmatism_map.jpg"

  proj="-JS0/-90/10c"
  rgn="-R0/360/-90/-60"

  gmt gmtset COLOR_NAN gray PS_MEDIA a0 MAP_FRAME_TYPE plain
  LABgrid=${LAB_mean}.grd
  gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
  gmt makecpt -T0/20/1 -Ccopper -I > age.cpt
  gmt grdimage $LABgrid $proj $rgn -Clitho.cpt -B30f10g30 -K -E600 -X12.5c -Y25c > $ps
  gmt pscoast $proj $rgn -Dh -W0.5p -A50000/0/0 -O -K >> $ps
  awk '{if ($4<=10) print $1, $2, $3}' ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean.xyaszs |\
    gmt psxy $proj $rgn -Cage.cpt -St0.2c -O -K >> $ps
  #awk -F, '{if(NR>1 && $5<65) print ($3 + $4)/2, ($1 + $2)/2, $5}' ${fold_data_output}/ANTARCTICA_rift_filtered.csv | gmt psxy $proj $rgn -Cage.cpt -Sc0.1c -O -K >> $ps
  #awk -F, '{if(NR>1 && $5<65) print $2, $1, $3}' ${fold_data_output}/BALL_filtered.csv | gmt psxy $proj $rgn -Cage.cpt -Sd0.1c -O -K >> $ps
  gmt psscale -D+5.125c/-1.5c+w10c/0.2c+jMC+h+e -Clitho.cpt -B50f10+l"z@-LAB@- (km)" -O -K >> $ps
  gmt psscale -D+5.125c/-3.5c+w10c/0.2c+jMC+h+e -Cage.cpt -B5f1+l"Age (Ma)" -O >> $ps

  gmt psconvert $ps -Tj -E600 -A0.1c/0.1c -P -Z
}

plot_LAB_magmatism_relationship(){

  gmt makecpt -T0/20/1 -Ccopper -I > age.cpt

  ps="./plots/LAB1200_${model}_magmatism_age.ps"
  jpg="./plots/LAB1200_${model}_magmatism_age.jpg"

  rgn="-R-0.5/23.5/0/200"
  scale="-JX21.5c/-10c"

  gmt psbasemap $rgn $scale -Bpx5f1+l"Age (Ma)" -Bpy20f10+l"z@-LAB@- (km)" -BWSne -K -Y25c > $ps
  awk '{if($4<=10) print $3, $5, $3, $4, $6}' ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean.xyaszs | gmt psxy $rgn $scale -E -St10p -O -Cage.cpt >> $ps
  gmt psconvert $ps -Tj -E600 -A0.1c -P -Z
}

plot_LAB_magmatism_correlation_histogram(){

  ps="./plots/LAB1200_magmatism_age_histogram.ps"
  jpg="./plots/LAB1200_magmatism_age_histogram.jpg"

  rgn="-R-0.3/1.0/0/5"
  scale="-JX10c/10c"

  #-Bpy0.5f0.25+l"Probability Density"
  gmt psbasemap $rgn $scale -Bpx0.2f0.1+l"@~\162@~" -Bpy0.5f0.25+l"f(@~\162@~)" -BWSne -K -Y25c > $ps
  awk '{print $3}' Antarctica_magmatism_LAB_CAM2016_spearman.spsp | gmt pshistogram -Z1 $rgn $scale -W0.01 -Gred -O -K -t50 >> $ps
  awk '{print $3}' Antarctica_magmatism_LAB_ANT-20_spearman.spsp | gmt pshistogram -Z1 $rgn $scale -W0.01 -Ggreen4 -O -K -t50 >> $ps
  awk '{print $3}' Antarctica_magmatism_LAB_SL2013_spearman.spsp | gmt pshistogram -Z1 $rgn $scale -W0.01 -Gblue -O -K -t50 >> $ps
  printf "0.296, 0.0\n0.296,5.0" | gmt psxy $rgn $scale -W1.0p,dashed -O >> $ps
  gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_scripts=${fold_base}/scripts
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')

fold_input_magmatism=$(awk '$1 ~ /^input_magmatism/' config.ini | awk '{print $3}')
fold_input_magmatism=${fold_data_input}/${fold_input_magmatism}
fold_input_LAB=$(awk '$1 ~ /^input_LAB_models/' config.ini | awk '{print $3}')
fold_input_LAB=${fold_data_input}/${fold_input_LAB}
fold_input_conductive_isotherms=$(awk '$1 ~ /^input_conductive_isotherms/' config.ini | awk '{print $3}')
fold_input_conductive_isotherms=${fold_data_input}/${fold_input_magmatism}

fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')

mkdir -p ${fold_plot_output}
mkdir -p ${fold_data_output}

model="ANT-20"
LAB_mean="${fold_input_LAB}/${model}/mean"
LAB_std="${fold_input_LAB}/${model}/std"
process_magmatism_data
#plot_LAB_magmatism_map
#plot_LAB_magmatism_relationship
model="CAM2016"
LAB_mean="${fold_input_LAB}/${model}/mean"
LAB_std="${fold_input_LAB}/${model}/std"
process_magmatism_data
#plot_LAB_magmatism_map
#plot_LAB_magmatism_relationship
model="SL2013"
LAB_mean="${fold_input_LAB}/${model}/mean"
LAB_std="${fold_input_LAB}/${model}/std"
process_magmatism_data
#plot_LAB_magmatism_map
#plot_LAB_magmatism_relationship
#plot_LAB_magmatism_correlation_histogram

rm -f gmt.conf gmt.history *.cpt junk*
