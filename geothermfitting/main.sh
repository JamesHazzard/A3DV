#!/bin/bash

copy_config(){

config="./config/jameshome.ini" # choose a config file
cp ${config} scripts/config.ini # copy config file to scripts
echo "date = ${date}" >> scripts/config.ini

}

make_output_dir(){

  python3 make_output_dir.py

}

remove_empty_dirs(){

  for entry in ${fold_logs}/*; do
    date_entry=$(echo "${entry##*/}")
    n_lines=$(cat $entry | wc -l)
    empty=$(echo $n_lines empty | awk '{if ($1 > 0) {print $1} else {print $2}}')
    if [ $empty == "empty" ]; then
      tmp_fold=${date_entry%%.*}
      rm -rf ${fold_data_output_archive}/${tmp_fold}
      rm -f ${fold_logs}/${date_entry}
      echo "deleting empty data output ${tmp_fold} and associated log file ${date_entry}"
    fi
  done

}

run_LAB1200_distribution(){

  make_output_dir

    for ((i=0; i<=999; i++)); do
        echo "running anelasticity model ${i}..."
        python3 make_LAB1200.py -v distribution -i ${i} # run script 
        echo "python3 make_LAB1200.py -v distribution -i ${i}" >> ${file_log} # update log with command used to run script
    done

    grid_type=LAB1200
    calculate_model_summary

}

run_LAB1200_MAP(){

  make_output_dir

  python3 make_LAB1200.py -v MAP # run script
  echo "python3 make_LAB1200.py -v MAP" >> ${file_log} # update log with command used to run script 

}

run_geotherms_MAP(){

  make_output_dir

  python3 make_geotherms.py -v MAP # run script
  #echo "python3 make_geotherms.py -v MAP" >> ${file_log} # update log with command used to run script

}

plot_geotherms_MAP(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  mkdir -p ${fold_plot_output}
  fold_plot_output=${fold_plot_output}/Tp_${potential_temperature}_sol50_${solidus_50km}
  mkdir -p ${fold_plot_output}
  fold_plot_output=${fold_plot_output}/geotherms
  mkdir -p ${fold_plot_output}

  f_locs=${fold_data_input}/${input_locs}/sample_locs.txt
  n_locs=8

  for ((i=1; i<=$n_locs; i++)); do

    lon=$(awk '{if(NR=='${i}') print $1}' ${f_locs})
    lat=$(awk '{if(NR=='${i}') print $2}' ${f_locs})
    f_geo=${fold_data_output}/Tp_${potential_temperature}_sol50_${solidus_50km}/latest/geotherms/MAP/lon_${lon}_lat_${lat}_geoth

    f_name=${fold_plot_output}/lon_${lon}_lat_${lat}_MAP_geotherm
    ps=${f_name}.ps
    jpg=${f_name}.jpg

    rgn="-R-10/1810/-10/410"
    scale="-JX10c/-21.5c"

    gmt psbasemap $rgn $scale -B0 -K -X25c -Y25c > $ps
    awk '{print $2, $1}' ${f_geo}_1.txt |\
      gmt psxy $rgn $scale -Sc0.18 -Ggrey -Wblack -O -K >> $ps # raw data
    #awk '{print $2, $1}' ${f_geo}_2.txt |\
    #  gmt psxy $rgn $scale -Sc0.18 -Gblack -Wblack -O -K >> $ps # crust removed
    #awk '{print $2, $1}' ${f_geo}_3.txt |\
    #  gmt psxy $rgn $scale -Sc0.18 -Gblue -Wblack -O -K >> $ps # spurious shallow dT/dz removed
    awk '{print $2, $1}' ${f_geo}_4.txt |\
      gmt psxy $rgn $scale -Sc0.18 -Gblack -Wblack -O -K >> $ps # basement reference added
    #awk '{print $2, $1}' ${f_geo}_5.txt |\
    #  gmt psxy $rgn $scale -Sc0.18 -Gblue -Wblack -O -K >> $ps # interpolated
    awk '{print $2, $1}' ${f_geo}_6.txt |\
      gmt psxy $rgn $scale -W2p,red -O -K >> $ps # fitted
    #awk '{print $2, $1}' "./geotherms/output_geotherms/geoth_${lon}_${lat}.zT" > junk
    #gmt psxy $rgn $scale junk -Sc0.1 -Gblue -Wblue -O -K >> $ps
    #awk '{print $2, $1}' "./geotherms/stitched_geotherms/geoth_${lon}_${lat}.zT" > junk
    #gmt psxy $rgn $scale junk -W4p,black -O -K -t25 >> $ps

    gmt psbasemap $rgn $scale -Bpx500f100+l"T (@~\260@~C)" -Bpy50f10+l"Depth (km)" -BWSne -O >> $ps
    gmt psconvert $ps -Tj -E600 -A0.1c/0.1c -P -Z

  done

}

plot_crustal_grids(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  mkdir -p ${fold_plot_output}
  fold_plot_output=${fold_plot_output}/Tp_${potential_temperature}_sol50_${solidus_50km}
  mkdir -p ${fold_plot_output}
  fold_plot_output=${fold_plot_output}/crustal_grids
  mkdir -p ${fold_plot_output}

  f_name="${fold_plot_output}/ANT-20_basement_depth"
  ps=${f_name}.ps
  jpg=${f_name}.jpg

  proj="-JS0/-90/10c"
  rgn="-R0/360/-90/-60"

  gmt gmtset COLOR_NAN gray PS_MEDIA a0 MAP_FRAME_TYPE plain
  gmt makecpt -T-2/2/0.25 -Cpolar -D -I > litho.cpt
  gmt grdimage ${fold_data_input}/crustal_grids/ant_cont_region_sampled_modified_crust1_basement_0.5d.grd $proj $rgn -Clitho.cpt -Bx15g15 -By1g1 -K -E300 -X12.5c -Y25c > $ps
  gmt pscoast $proj $rgn -Dh -W0.5p -A50000/0/0 -O -K >> $ps
  gmt psscale -D+5.125c/-1.5c+w10c/0.2c+jMC+h+e -Clitho.cpt -B2f1+l"Basement depth (km)" -O >> $ps

  gmt psconvert $ps -Tj -E300 -A0.5c/0.5c -P -Z

}

calculate_model_summary(){

  vsmodel=ANT-20

  fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
  fold_previous_date="2023-01-13_11-47-37" # specify a previous run date to calculate model summary based on
  potential_temperature=$(awk '$1 ~ /^potential_temperature/' config.ini | awk '{print $3}')  # ensure vars in config match previous run date
  solidus_50km=$(awk '$1 ~ /^solidus_50km/' config.ini | awk '{print $3}')
  year=$(echo $fold_previous_date | tr "_" "\t" | awk '{print $1}' | tr "-" "\t" | awk '{print $1}')
  month=$(echo $fold_previous_date | tr "_" "\t" | awk '{print $1}' | tr "-" "\t" | awk '{print $2}')
  day=$(echo $fold_previous_date | tr "_" "\t" | awk '{print $1}' | tr "-" "\t" | awk '{print $3}')

  if [ ${grid_type} == baseTBL_const_Tp ]; then

    grid_type_in=HF_const_Tp
    fold_grid_in=${fold_data_output_archive}/RUN_${fold_previous_date}/${grid_type_in}/distribution
    fold_grid_out=${fold_data_output_archive}/RUN_${fold_previous_date}/${grid_type}/distribution
    mkdir -p ${fold_grid_out}

  else

    grid_type_in=${grid_type}
    fold_grid=${fold_data_output_archive}/RUN_${fold_previous_date}/${grid_type}/distribution
    fold_grid_in=${fold_grid}
    fold_grid_out=${fold_grid}

  fi

  fold_grid_ind_in=${fold_grid_in}/individual
  fold_grid_ind_out=${fold_grid_out}/individual
  fold_grid_summary=${fold_grid_out}/summary
  mkdir -p ${fold_grid_ind_out}
  mkdir -p ${fold_grid_summary}

  rgn="-R0/360/-90/-60"
  inc="-I0.5d/0.5d"

  n_files=$(ls ${fold_grid_ind_in}/*.txt | wc -l)
  idx_max=$(echo $n_files | awk '{print $1-1}')

  idx=0
  echo mean $idx
  path_in=${fold_grid_ind_in}/${vsmodel}_${grid_type_in}_${idx}
  path_out=${fold_grid_ind_out}/${vsmodel}_${grid_type}_${idx}

  if [ ${grid_type} == HF_const_Tp ] || [ ${grid_type} == HF ]; then
    echo "picking ${grid_type} output from .txt files and converting to grid.."
    conversion_call=$(awk '{print $1, $2, $12}' ${path_in}.txt | gmt xyz2grd $rgn $inc -G${path_out}.grd)
    $conversion_call
  elif [ ${grid_type} == baseTBL_const_Tp ]; then
    echo "picking baseTBL output from .txt files and converting to grid.."
    awk '{print $1, $2, $10}' ${path_in}.txt > ${path_out}.txt
    conversion_call=$(awk '{print $1, $2, $10}' ${path_in}.txt | gmt xyz2grd $rgn $inc -G${path_out}.grd)
    $conversion_call
  else
    echo "converting .txt file to grid.."
    conversion_call=$(gmt xyz2grd $rgn $inc ${path_in}.txt -G${path_out}.grd)
    $conversion_call
  fi

  gmt grdmath ${path_out}.grd 0 DENAN = tmp_${grid_type}.grd
  gmt grdmath tmp_${grid_type}.grd 0 NEQ = tmp_data_count_${grid_type}.grd

  cp tmp_${grid_type}.grd sum_${grid_type}.grd
  cp tmp_data_count_${grid_type}.grd sum_data_count_${grid_type}.grd

  for ((idx=1; idx<=$idx_max; idx++)); do

    echo mean $idx
    path_in=${fold_grid_ind_in}/${vsmodel}_${grid_type_in}_${idx}
    path_out=${fold_grid_ind_out}/${vsmodel}_${grid_type}_${idx}

    if [ ${grid_type} == HF_const_Tp ] || [ ${grid_type} == HF ]; then
      echo "picking HF output from .txt files and converting to grid.."
      conversion_call=$(awk '{print $1, $2, $12}' ${path_in}.txt | gmt xyz2grd $rgn $inc -G${path_out}.grd)
      $conversion_call
    elif [ ${grid_type} == baseTBL_const_Tp ]; then
      echo "picking baseTBL output from .txt files and converting to grid.."
      awk '{print $1, $2, $10}' ${path_in}.txt > ${path_out}.txt
      conversion_call=$(awk '{print $1, $2, $10}' ${path_in}.txt | gmt xyz2grd $rgn $inc -G${path_out}.grd)
      $conversion_call
    else
      echo "converting .txt file to grid.."
      conversion_call=$(gmt xyz2grd $rgn $inc ${path_in}.txt -G${path_out}.grd)
      $conversion_call
    fi
    gmt grdmath ${path_out}.grd 0 DENAN = tmp_${grid_type}.grd
    gmt grdmath tmp_${grid_type}.grd 0 NEQ = tmp_data_count_${grid_type}.grd

    gmt grdmath tmp_${grid_type}.grd sum_${grid_type}.grd ADD = sum_${grid_type}.grd
    gmt grdmath tmp_data_count_${grid_type}.grd sum_data_count_${grid_type}.grd ADD = sum_data_count_${grid_type}.grd

  done

  gmt grdmath sum_${grid_type}.grd sum_data_count_${grid_type}.grd DIV = mean_${grid_type}.grd


  idx=0
  echo std $idx
  path_in=${fold_grid_ind_in}/${vsmodel}_${grid_type_in}_${idx}
  path_out=${fold_grid_ind_out}/${vsmodel}_${grid_type}_${idx}

  gmt grdmath ${path_out}.grd 0 DENAN = tmp_${grid_type}.grd
  gmt grdmath tmp_${grid_type}.grd mean_${grid_type}.grd SUB SQR = dist_${grid_type}.grd

  for ((idx=1; idx<=$idx_max; idx++)); do

    echo std $idx
    path_in=${fold_grid_ind_in}/${vsmodel}_${grid_type_in}_${idx}
    path_out=${fold_grid_ind_out}/${vsmodel}_${grid_type}_${idx}
    gmt grdmath ${path_out}.grd 0 DENAN = tmp_${grid_type}.grd
    gmt grdmath tmp_${grid_type}.grd mean_${grid_type}.grd SUB SQR = tmp_dist_${grid_type}.grd
    gmt grdmath tmp_dist_${grid_type}.grd dist_${grid_type}.grd ADD = dist_${grid_type}.grd

  done

  gmt grdmath dist_${grid_type}.grd sum_data_count_${grid_type}.grd DIV SQRT = std_${grid_type}.grd

  cp sum_data_count_${grid_type}.grd ${fold_grid_summary}/${grid_type}_data_count.grd
  cp mean_${grid_type}.grd ${fold_grid_summary}/${grid_type}_mean.grd
  cp std_${grid_type}.grd ${fold_grid_summary}/${grid_type}_std.grd
  gmt grdfilter ${fold_grid_summary}/${grid_type}_mean.grd -D4 -Fg400 -G${fold_grid_summary}/${grid_type}_mean_filtered.grd
  gmt grd2xyz ${fold_grid_summary}/${grid_type}_mean_filtered.grd > ${fold_grid_summary}/${grid_type}_mean_filtered.xyz
  gmt grdfilter ${fold_grid_summary}/${grid_type}_std.grd -D4 -Fg400 -G${fold_grid_summary}/${grid_type}_std_filtered.grd

  cp -r ${fold_grid_out} ${fold_data_output_latest}/${grid_type}

  rm tmp*.grd sum*.grd mean*.grd std*.grd dist*.grd

}

filter_model_summary(){

  #fold_previous_date="2022-12-21_15-49-23" # specify a previous run date to filter model summary on
  fold_previous_date="2022-12-09_17-13-56"
  fold_grid=${fold_data_output_archive}/RUN_${fold_previous_date}/${grid_type}/distribution
  fold_grid_summary=${fold_grid}/summary
  gmt grdfilter ${fold_grid_summary}/${grid_type}_mean.grd -D4 -Fg400 -G${fold_grid_summary}/${grid_type}_mean_filtered.grd
  gmt grd2xyz ${fold_grid_summary}/${grid_type}_mean_filtered.grd > ${fold_grid_summary}/${grid_type}_mean_filtered.xyz
  gmt grdfilter ${fold_grid_summary}/${grid_type}_std.grd -D4 -Fg400 -G${fold_grid_summary}/${grid_type}_std_filtered.grd
  cp -r ${fold_grid} ${fold_data_output_latest}/${grid_type}

}

get_hpc_data(){

  # write your script here to pull data from your hpc (high-performance computing, aka cluster, aka supercomputer) server into your data output space
  # this is left incomplete here, in order to protect my credentials, and because file transfer protocols vary
  
  make_output_dir
  grid_type="HF"
  echo "this data was generated on the HPC and involves the script make_${grid_type}.py" >> ${file_log}
  mkdir -p ${fold_data_output_archive}/RUN_${fold_date}/${grid_type}
  mkdir -p ${fold_data_output_archive}/RUN_${fold_date}/${grid_type}/distribution
  mkdir -p ${fold_data_output_archive}/RUN_${fold_date}/${grid_type}/distribution/individual
  mkdir -p ${fold_data_output_archive}/RUN_${fold_date}/${grid_type}/distribution/summary

}

make_latest_output_dirs(){

  for grid_type in LAB1200 HF_const_Tp HF baseTBL_const_Tp geotherms; do

    mkdir -p ${fold_data_output_latest}
    mkdir -p ${fold_data_output_latest}/${grid_type}
    mkdir -p ${fold_data_output_latest}/${grid_type}/distribution
    mkdir -p ${fold_data_output_latest}/${grid_type}/distribution/individual
    mkdir -p ${fold_data_output_latest}/${grid_type}/distribution/summary
    mkdir -p ${fold_data_output_latest}/${grid_type}/MAP

  done

  mkdir -p ${fold_data_output_latest}/LAB1200_to_baseTBL

}

date=$(date '+%Y-%m-%d_%H-%M-%S')
vsmodel="ANT-20"

copy_config
cd scripts  # change working directory to location of scripts
make all
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
input_locs=$(awk '$1 ~ /^input_locs/' config.ini | awk '{print $3}')
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}') # create log file
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')
fold_date=$(awk '$1 ~ /^date/' config.ini | awk '{print $3}')
potential_temperature=$(awk '$1 ~ /^potential_temperature/' config.ini | awk '{print $3}')
solidus_50km=$(awk '$1 ~ /^solidus_50km/' config.ini | awk '{print $3}')
fold_data_output_archive=${fold_data_output}/Tp_${potential_temperature}_sol50_${solidus_50km}/archive
fold_data_output_latest=${fold_data_output}/Tp_${potential_temperature}_sol50_${solidus_50km}/latest
fold_logs=${fold_data_output_archive}/logs

file_log=${fold_logs}/RUN_${fold_date}.log
echo -n >> ${file_log}

#get_hpc_data
#make_latest_output_dirs
#grid_type="HF"
#calculate_model_summary
#run_geotherms_MAP
plot_geotherms_MAP
#plot_crustal_grids

rm -f gmt.* input.dat output.dat geoth.out lacdevs.dat
remove_empty_dirs