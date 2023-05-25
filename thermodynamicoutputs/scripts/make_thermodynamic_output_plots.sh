#!/bin/bash

plot_thermodynamic_distance_histograms(){

  # Set GMT plotting parameters   
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  compare_Tp=$(awk '$1 ~ /^compare_potential_temperature/' config.ini | awk '{print $3}')
  compare_sol50=$(awk '$1 ~ /^compare_solidus_50km/' config.ini | awk '{print $3}')
  compare_date=$(awk '$1 ~ /^compare_date/' config.ini | awk '{print $3}')

  fold_actual_plot_output=${fold_plot_output}/comparison
  mkdir -p ${fold_actual_plot_output}

  # make a new output dir if this combination of dates has not been looked at before
  # if it has, still make the plot but output it to the pre-exisiting dir (overwrite)
  output_digit=$(find ${fold_actual_plot_output} -mindepth 1 -maxdepth 1 -type d | wc -l | awk '{print $1+1}' | awk '{printf "%.2d", $1}')
  outputs=$(find ${fold_actual_plot_output} -mindepth 1 -maxdepth 1 -type d | sort)

  for output in $outputs; do

    find_date=$(grep -Rc "${date}" ${output}/comparison.log)
    find_compare_date=$(grep -Rc "${compare_date}" ${output}/comparison.log)
    find_both_dates=$(echo ${find_date} ${find_compare_date} | awk '{print $1*$2}')

    if [ $find_both_dates -eq 1 ]; then

      fold_actual_plot_output=${output}

    else

      fold_actual_plot_output=${fold_actual_plot_output}/comparison_${output_digit}
      mkdir -p ${fold_actual_plot_output}
      f_log=${fold_actual_plot_output}/comparison.log
      echo -n > $f_log
      echo date_1 ${date} >> $f_log
      echo Tp_1 ${Tp} >> $f_log
      echo sol50_1 ${sol50} >> $f_log
      echo date_2 ${compare_date} >> $f_log
      echo Tp_2 ${compare_Tp} >> $f_log
      echo sol50_2 ${compare_sol50} >> $f_log

    fi

  done

  # now we have an output dir, set up the plotting
  f_name="thermo_outputs_distance_histogram_panel"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  rgn="-R0/1/0/1"
  scale="-JX10c/10c"
  gmt psbasemap $rgn $scale -Y30c -B0 -K > $ps

  # global temperature comparison
  grid_type="temperature"
  full_grid_type="potential_temperature"
  fold_grid=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  fold_grid_compare=${fold_data_output}/${compare_date}/Tp_${compare_Tp}_sol50_${compare_sol50}/${grid_type}_grids/${summary_type}

  for depth in 75 150 250 350; do

    x_label="@~\161@~ / 1000 (@~\260@~C, ${depth} km)"
    grid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean
    grid_compare=${fold_grid_compare}/${full_grid_type}_depth_${depth}_Tp_${compare_Tp}_sol50_${compare_sol50}_mean
    grid_distance=${fold_actual_plot_output}/${full_grid_type}_depth_${depth}_mean_distance
    x_major_tick=0.1
    x_minor_tick=0.025
    bin_width=0.025
    echo "working on depth slice ${depth} km!"
    echo "checking for .xyz files.."
    if [[ ! -e ${grid_distance}.xyz ]]; then
      gmt grdmath ${grid}.grd ${grid_compare}.grd SUB = ${grid_distance}.grd
      gmt grd2xyz ${grid_distance}.grd > ${grid_distance}.xyz
    fi
    echo "calculating distribution ranges.."
    #awk '{print $3}' ${grid}.xyz | gmt pshistogram -Z1 ${scale} -W${bin_width} -I > junk1
    #x_min=$(awk '{print $1}' junk1)
    #x_max=$(awk '{print $2}' junk1)
    #y_max=$(awk '{print 1.05*$4}' junk1)
    #awk '{print $3}' ${grid_compare}.xyz | gmt pshistogram -Z1 ${scale} -W${bin_width} -I > junk2
    #x_min=$(awk '{print $1, '${x_min}'}' junk2 | awk '{if ($1 < $2) print $1; else print $2}')
    #x_min=0
    #x_max=$(awk '{print $2, '${x_max}'}' junk2 | awk '{if ($1 > $2) print $1; else print $2}')
    #x_max=1.8
    #y_max=$(awk '{print 1.05*$4, '${y_max}'}' junk2 | awk '{if ($1 > $2) print $1; else print $2}')
    x_min=-0.25
    x_max=0.25
    y_max=50
    echo $x_min $x_max $y_max
    rgn="-R${x_min}/${x_max}/0/${y_max}"
    if [[ ${depth} -eq 75 ]]; then
      move=""
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -Bpy10f5+l"Frequency (%)" -BWSne $move -O -K >> $ps
    else
      move="-X11c"
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -BwSne $move -O -K >> $ps
    fi
    echo "plotting.."
    awk '{print $3}' ${grid_distance}.xyz | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gblack -O -K -t10 >> $ps
    echo "done"

  done

  gmt psbasemap $rgn $scale -B0 -X-33c -Y-12c -O -K >> $ps

  # global viscosity comparison
  grid_type="viscosity"
  full_grid_type="viscosity"
  fold_grid=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  fold_grid_compare=${fold_data_output}/${compare_date}/Tp_${compare_Tp}_sol50_${compare_sol50}/${grid_type}_grids/${summary_type}

  for depth in 75 150 250 350; do

    x_label="log@-10@-@~\155@~@-@~\150@~@- (Pa s, ${depth} km)"
    grid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean
    grid_compare=${fold_grid_compare}/${full_grid_type}_depth_${depth}_Tp_${compare_Tp}_sol50_${compare_sol50}_mean
    grid_distance=${fold_actual_plot_output}/${full_grid_type}_depth_${depth}_mean_distance
    x_major_tick=1
    x_minor_tick=0.25
    bin_width=0.25
    echo "working on depth slice ${depth} km!"
    echo "checking for .xyz files.."
    rerun=$(echo 0)

    if [[ ! -e ${grid_distance}.xyz || ! -e ${grid_distance}_no_lithosphere.xyz || rerun -eq 1 ]]; then

      echo calculating distance between grids..
      gmt grdmath ${grid}.grd 22.5 LT = junk1.grd # 1 if less than 22.5, 0 otherwise
      gmt grdmath 1 junk1.grd DIV = junk2.grd # 1 if less than 22.5, NaN otherwise
      gmt grdmath ${grid}.grd junk2.grd MUL = junk3.grd # ${value} if less than 22.5, NaN otherwise
      gmt grdmath ${grid_compare}.grd 22.5 LT = junk4.grd # 1 if less than 22.5, 0 otherwise
      gmt grdmath 1 junk4.grd DIV = junk5.grd # 1 if less than 22.5, NaN otherwise
      gmt grdmath ${grid_compare}.grd junk5.grd MUL = junk6.grd # ${value} if less than 22.5, NaN otherwise
      gmt grdmath junk3.grd junk6.grd SUB = ${grid_distance}_no_lithosphere.grd
      gmt grdmath ${grid}.grd ${grid_compare}.grd SUB = ${grid_distance}.grd

      gmt grd2xyz ${grid_distance}.grd > ${grid_distance}.xyz
      gmt grd2xyz ${grid_distance}_no_lithosphere.grd > ${grid_distance}_no_lithosphere.xyz

    fi
    echo "calculating distribution ranges.."
    x_min=-5
    x_max=5
    y_max=50
    echo $x_min $x_max $y_max
    rgn="-R${x_min}/${x_max}/0/${y_max}"
    if [[ ${depth} -eq 75 ]]; then
      move=""
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -Bpy10f5+l"Frequency (%)" -BWSne $move -O -K >> $ps
    else
      move="-X11c"
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -BwSne $move -O -K >> $ps
    fi
    echo "plotting.."
    awk '{print $3}' ${grid_distance}.xyz | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gblack -O -K -t10 >> $ps
    echo "done"

  done

  gmt psbasemap $rgn $scale -B0 -X-33c -Y-12c -O -K >> $ps
  
  # global viscosity comparison
  grid_type="viscosity"
  full_grid_type="viscosity"
  fold_grid=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  fold_grid_compare=${fold_data_output}/${compare_date}/Tp_${compare_Tp}_sol50_${compare_sol50}/${grid_type}_grids/${summary_type}

  for depth in 75 150 250 350; do

    x_label="log@-10@-@~\155@~@-@~\150@~@- (Pa s, ${depth} km)"
    grid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean
    grid_compare=${fold_grid_compare}/${full_grid_type}_depth_${depth}_Tp_${compare_Tp}_sol50_${compare_sol50}_mean
    grid_distance=${fold_actual_plot_output}/${full_grid_type}_depth_${depth}_mean_distance
    x_major_tick=0.2
    x_minor_tick=0.1
    bin_width=0.1
    echo "working on depth slice ${depth} km!"
    echo "checking for .xyz files.."
    rerun=$(echo 0)

    if [[ ! -e ${grid_distance}.xyz || ! -e ${grid_distance}_no_lithosphere.xyz || rerun -eq 1 ]]; then

      echo calculating distance between grids..
      gmt grdmath ${grid}.grd 22.5 LT = junk1.grd # 1 if less than 22.5, 0 otherwise
      gmt grdmath 1 junk1.grd DIV = junk2.grd # 1 if less than 22.5, NaN otherwise
      gmt grdmath ${grid}.grd junk2.grd MUL = junk3.grd # ${value} if less than 22.5, NaN otherwise
      gmt grdmath ${grid_compare}.grd 22.5 LT = junk4.grd # 1 if less than 22.5, 0 otherwise
      gmt grdmath 1 junk4.grd DIV = junk5.grd # 1 if less than 22.5, NaN otherwise
      gmt grdmath ${grid_compare}.grd junk5.grd MUL = junk6.grd # ${value} if less than 22.5, NaN otherwise
      gmt grdmath junk3.grd junk6.grd SUB = ${grid_distance}_no_lithosphere.grd
      gmt grdmath ${grid}.grd ${grid_compare}.grd SUB = ${grid_distance}.grd

      gmt grd2xyz ${grid_distance}.grd > ${grid_distance}.xyz
      gmt grd2xyz ${grid_distance}_no_lithosphere.grd > ${grid_distance}_no_lithosphere.xyz
      grid_distance=${grid_distance}_no_lithosphere

    fi
    echo "calculating distribution ranges.."
    x_min=0
    x_max=1
    y_max=50
    echo $x_min $x_max $y_max
    rgn="-R${x_min}/${x_max}/0/${y_max}"
    if [[ ${depth} -eq 75 ]]; then
      move=""
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -Bpy10f5+l"Frequency (%)" -BWSne $move -O -K >> $ps
    else
      move="-X11c"
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -BwSne $move -O -K >> $ps
    fi
    echo "plotting.."
    awk '{print $3}' ${grid_distance}.xyz | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gblack -O -K -t10 >> $ps
    echo "done"

  done

  gmt psbasemap $rgn $scale -B0 -O >> $ps
  gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

plot_thermodynamic_histograms(){

  # Set GMT plotting parameters   
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  compare_Tp=$(awk '$1 ~ /^compare_potential_temperature/' config.ini | awk '{print $3}')
  compare_sol50=$(awk '$1 ~ /^compare_solidus_50km/' config.ini | awk '{print $3}')
  compare_date=$(awk '$1 ~ /^compare_date/' config.ini | awk '{print $3}')

  fold_actual_plot_output=${fold_plot_output}/comparison
  mkdir -p ${fold_actual_plot_output}

  # make a new output dir if this combination of dates has not been looked at before
  # if it has, still make the plot but output it to the pre-exisiting dir (overwrite)
  output_digit=$(find ${fold_actual_plot_output} -mindepth 1 -maxdepth 1 -type d | wc -l | awk '{print $1+1}' | awk '{printf "%.2d", $1}')
  outputs=$(find ${fold_actual_plot_output} -mindepth 1 -maxdepth 1 -type d | sort)

  for output in $outputs; do

    find_date=$(grep -Rc "${date}" ${output}/comparison.log)
    find_compare_date=$(grep -Rc "${compare_date}" ${output}/comparison.log)
    find_both_dates=$(echo ${find_date} ${find_compare_date} | awk '{print $1*$2}')

    if [ $find_both_dates -eq 1 ]; then

      fold_actual_plot_output=${output}

    else

      fold_actual_plot_output=${fold_actual_plot_output}/comparison_${output_digit}
      mkdir -p ${fold_actual_plot_output}
      f_log=${fold_actual_plot_output}/comparison.log
      echo -n > $f_log
      echo date_1 ${date} >> $f_log
      echo Tp_1 ${Tp} >> $f_log
      echo sol50_1 ${sol50} >> $f_log
      echo date_2 ${compare_date} >> $f_log
      echo Tp_2 ${compare_Tp} >> $f_log
      echo sol50_2 ${compare_sol50} >> $f_log

    fi

  done

  # now we have an output dir, set up the plotting
  f_name="thermo_outputs_histogram_panel"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  rgn="-R0/1/0/1"
  scale="-JX10c/10c"
  gmt psbasemap $rgn $scale -Y30c -B0 -K > $ps

  # global temperature comparison
  grid_type="temperature"
  full_grid_type="potential_temperature"
  fold_grid=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  fold_grid_compare=${fold_data_output}/${compare_date}/Tp_${compare_Tp}_sol50_${compare_sol50}/${grid_type}_grids/${summary_type}

  for depth in 75 150 250 350; do

    x_label="@~\161@~ / 1000 (@~\260@~C, ${depth} km)"
    grid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean
    grid_compare=${fold_grid_compare}/${full_grid_type}_depth_${depth}_Tp_${compare_Tp}_sol50_${compare_sol50}_mean
    x_major_tick=0.2
    x_minor_tick=0.1
    bin_width=0.1
    echo "working on depth slice ${depth} km!"
    echo "checking for .xyz files.."
    if [[ ! -e ${grid}.xyz || ! -e ${grid_compare}.xyz ]]; then
      gmt grd2xyz ${grid}.grd > ${grid}.xyz
      gmt grd2xyz ${grid_compare}.grd > ${grid_compare}.xyz
    fi
    echo "calculating distribution ranges.."
    #awk '{print $3}' ${grid}.xyz | gmt pshistogram -Z1 ${scale} -W${bin_width} -I > junk1
    #x_min=$(awk '{print $1}' junk1)
    #x_max=$(awk '{print $2}' junk1)
    #y_max=$(awk '{print 1.05*$4}' junk1)
    #awk '{print $3}' ${grid_compare}.xyz | gmt pshistogram -Z1 ${scale} -W${bin_width} -I > junk2
    #x_min=$(awk '{print $1, '${x_min}'}' junk2 | awk '{if ($1 < $2) print $1; else print $2}')
    #x_min=0
    #x_max=$(awk '{print $2, '${x_max}'}' junk2 | awk '{if ($1 > $2) print $1; else print $2}')
    #x_max=1.8
    #y_max=$(awk '{print 1.05*$4, '${y_max}'}' junk2 | awk '{if ($1 > $2) print $1; else print $2}')
    x_min=0
    x_max=1.8
    y_max=50
    echo $x_min $x_max $y_max
    rgn="-R${x_min}/${x_max}/0/${y_max}"
    if [[ ${depth} -eq 75 ]]; then
      move=""
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -Bpy10f5+l"Frequency (%)" -BWSne $move -O -K >> $ps
    else
      move="-X11c"
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -BwSne $move -O -K >> $ps
    fi
    echo "plotting.."
    awk '{print $3}' ${grid}.xyz | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gblue -O -K -t50 >> $ps
    awk '{print $3}' ${grid_compare}.xyz | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gred -O -K -t50 >> $ps
    echo "done"

  done

  gmt psbasemap $rgn $scale -B0 -X-33c -Y-12c -O -K >> $ps

  # global viscosity comparison
  grid_type="viscosity"
  full_grid_type="viscosity"
  fold_grid=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  fold_grid_compare=${fold_data_output}/${compare_date}/Tp_${compare_Tp}_sol50_${compare_sol50}/${grid_type}_grids/${summary_type}

  for depth in 75 150 250 350; do

    x_label="log@-10@-@~\155@~@-@~\150@~@- (Pa s, ${depth} km)"
    grid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean
    grid_compare=${fold_grid_compare}/${full_grid_type}_depth_${depth}_Tp_${compare_Tp}_sol50_${compare_sol50}_mean
    x_major_tick=2
    x_minor_tick=0.5
    bin_width=0.5
    echo "working on depth slice ${depth} km!"
    echo "checking for .xyz files.."
    if [[ ! -e ${grid}.xyz || ! -e ${grid_compare}.xyz ]]; then
      gmt grd2xyz ${grid}.grd > ${grid}.xyz
      gmt grd2xyz ${grid_compare}.grd > ${grid_compare}.xyz
    fi
    echo "calculating distribution ranges.."
    #awk '{print $3}' ${grid}.xyz | gmt pshistogram -Z1 ${scale} -W${bin_width} -I > junk1
    #x_min=$(awk '{print $1}' junk1)
    #x_max=$(awk '{print $2}' junk1)
    #y_max=$(awk '{print 1.05*$4}' junk1)
    #awk '{print $3}' ${grid_compare}.xyz | gmt pshistogram -Z1 ${scale} -W${bin_width} -I > junk2
    #x_min=$(awk '{print $1, '${x_min}'}' junk2 | awk '{if ($1 < $2) print $1; else print $2}')
    #x_min=0
    #x_max=$(awk '{print $2, '${x_max}'}' junk2 | awk '{if ($1 > $2) print $1; else print $2}')
    #x_max=1.8
    #y_max=$(awk '{print 1.05*$4, '${y_max}'}' junk2 | awk '{if ($1 > $2) print $1; else print $2}')
    x_min=17
    x_max=30
    y_max=50
    echo $x_min $x_max $y_max
    rgn="-R${x_min}/${x_max}/0/${y_max}"
    if [[ ${depth} -eq 75 ]]; then
      move=""
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -Bpy10f5+l"Frequency (%)" -BWSne $move -O -K >> $ps
    else
      move="-X11c"
      gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -BwSne $move -O -K >> $ps
    fi
    echo "plotting.."
    awk '{print $3}' ${grid}.xyz | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gred -O -K -t50 >> $ps
    awk '{print $3}' ${grid_compare}.xyz | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gblue -O -K -t50 >> $ps
    echo "done"

  done

  
  gmt psbasemap $rgn $scale -B0 -O >> $ps
  gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

convert_temperature_potential_temperature_grids(){

  grid_type="temperature"
  full_grid_type="potential_temperature"
  fold_grid=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}

  rgn="-R0/360/-90/-60"
  grid_increment="-I0.05d/0.05d"

  alpha=0.00003
  Cp=1187
  g=9.81
  e=2.718281828459045

  for depth in 75 150 250 350; do

    Tp_m=$(echo $alpha $g $depth $Cp $e | awk '{print $5^(-($1*$2*$3*1000)/$4)}')
    echo $Tp_m
    Tgrid_mean=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
    Tgrid_std=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
    Tpgrid_mean=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
    Tpgrid_std=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd

    if [[ -e $Tpgrid_mean ]]; then
      echo "file exists, skipping.."
    else
      gmt grdmath $Tgrid_mean 273.15 ADD ${Tp_m} MUL 273.15 SUB 1000 DIV = $Tpgrid_mean
      gmt grdmath $Tgrid_std ${Tp_m} MUL = $Tpgrid_std
    fi

    echo depth slice $depth km done

  done

}

convert_temperature_potential_temperature_comparison_grids(){

  grid_type="temperature_comparison"
  full_grid_type="potential_temperature_comparison"
  fold_grid=${fold_actual_data_output}/temperature_comparison_grids/distribution/summary

  cd ${fold_grid}

  rgn="-R0/360/-90/-60"
  grid_increment="-I0.05d/0.05d"

  alpha=0.00003
  Cp=1187
  g=9.81
  e=2.718281828459045

  for depth in 75 150 250 350; do

    Tp_m=$(echo $alpha $g $depth $Cp $e | awk '{print $5^(-($1*$2*$3*1000)/$4)}')
    echo $Tp_m
    Tgrid_mean=${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
    Tgrid_std=${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
    Tpgrid_mean=${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
    Tpgrid_std=${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd

    if [[ -e $Tpgrid_mean ]]; then
      echo "file exists, skipping.."
    else
      gmt grdmath $Tgrid_mean 273.15 ADD ${Tp_m} MUL 273.15 SUB 1000 DIV = $Tpgrid_mean
      gmt grdmath $Tgrid_std ${Tp_m} MUL = $Tpgrid_std
    fi

    echo depth slice $depth km done

  done

  cd ${fold_scripts}

}

get_comparison_plot_output(){

    compare_Tp=$(awk '$1 ~ /^compare_potential_temperature/' config.ini | awk '{print $3}')
    compare_sol50=$(awk '$1 ~ /^compare_solidus_50km/' config.ini | awk '{print $3}')
    compare_date=$(awk '$1 ~ /^compare_date/' config.ini | awk '{print $3}')

    fold_actual_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')/comparison
    fold_actual_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')/comparison
    mkdir -p ${fold_actual_plot_output}

    # make a new output dir if this combination of dates has not been looked at before
    # if it has, still make the plot but output it to the pre-exisiting dir (overwrite)
    output_digit=$(find ${fold_actual_plot_output} -mindepth 1 -maxdepth 1 -type d | wc -l | awk '{print $1+1}' | awk '{printf "%.2d", $1}')
    outputs=$(find ${fold_actual_plot_output} -mindepth 1 -maxdepth 1 -type d | sort)

    if [[ $(echo ${outputs} | awk '{print NF}' | awk '{if($1<1) print 999}') -eq 999 ]]; then
        
        fold_actual_plot_output=${fold_actual_plot_output}/comparison_${output_digit}
        mkdir -p ${fold_actual_plot_output}
        f_log=${fold_actual_plot_output}/comparison.log
        echo -n > $f_log
        echo date_1 ${date} >> $f_log
        echo Tp_1 ${Tp} >> $f_log
        echo sol50_1 ${sol50} >> $f_log
        echo date_2 ${compare_date} >> $f_log
        echo Tp_2 ${compare_Tp} >> $f_log
        echo sol50_2 ${compare_sol50} >> $f_log

    else

        for output in $outputs; do

            find_date=$(grep -Rc "${date}" ${output}/comparison.log)
            find_compare_date=$(grep -Rc "${compare_date}" ${output}/comparison.log)
            find_both_dates=$(echo ${find_date} ${find_compare_date} | awk '{print $1*$2}')

            if [ $find_both_dates -eq 1 ]; then

                fold_actual_plot_output=${output}
                output_digit=${output: -2}

            else

                fold_actual_plot_output=${fold_actual_plot_output}/comparison_${output_digit}
                mkdir -p ${fold_actual_plot_output}
                f_log=${fold_actual_plot_output}/comparison.log
                echo -n > $f_log
                echo date_1 ${date} >> $f_log
                echo Tp_1 ${Tp} >> $f_log
                echo sol50_1 ${sol50} >> $f_log
                echo date_2 ${compare_date} >> $f_log
                echo Tp_2 ${compare_Tp} >> $f_log
                echo sol50_2 ${compare_sol50} >> $f_log

            fi

        done

    fi

    fold_actual_data_output=${fold_actual_data_output}/comparison_${output_digit}

}

plot_temperature_distance_quad_panel(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  grid_type="temperature_comparison"
  full_grid_type="potential_temperature_comparison"
  fold_grid=${fold_actual_data_output}/temperature_comparison_grids/distribution/summary

  f_name="${grid_type}_quad_panel_mean"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  gmt makecpt -Cpolar -D -T-0.3/0.3/0.03 > temp.cpt

  depth=75
  # Plot panel a, mean temperature at 75 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y32c -X20c > $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 75 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps

  depth=150
  # Plot panel b, mean temperature at 150 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx-4.5c/-1.1c+w7c/0.25c+e+h -Ctemp.cpt -B0.1f0.05+l"@~\155@~@-@~\104\161@~@- / 1000 (@~\260@~C)" -Al -O -K >> $ps
  echo "0.6 1.05 150 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps

  depth=250
  # Plot panel c, mean temperature at 250 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y-12.5c -X-12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 250 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 c" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps

  depth=350
  # Plot panel d, mean temperature at 350 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "0.6 1.05 350 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 d" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps

  gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
  rm *.cpt

  # Set plot name
  grid_type="temperature_comparison"
  full_grid_type="potential_temperature_comparison"
  fold_grid=${fold_actual_data_output}/temperature_comparison_grids/distribution/summary

  f_name="${grid_type}_quad_panel_std"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  gmt makecpt -Cpolar -D -T0/60/3 > temp.cpt

  depth=75
  # Plot panel a, std temperature at 75 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y32c -X10c > $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 75 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=150
  # Plot panel b, std temperature at 150 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx-4.5c/-1.1c+w7c/0.25c+e+h -Ctemp.cpt -B10f5+l"@~\163@~@-@~\104\161@~@- (@~\260@~C)" -Al -O -K >> $ps
  echo "0.6 1.05 150 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=250
  # Plot panel c, std temperature at 250 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y-12.5c -X-12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 250 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 c" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=350
  # Plot panel d, std temperature at 350 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "0.6 1.05 350 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 d" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps

  gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
  rm *.cpt

  cd ${fold_base}/scripts

}

plot_viscosity_distance_quad_panel(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  grid_type="viscosity_comparison"
  fold_grid=${fold_actual_data_output}/viscosity_comparison_grids/distribution/summary

  f_name="${grid_type}_quad_panel_mean"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  gmt makecpt -Cpolar -D -T-2/2/0.25 -I > visc.cpt

  depth=75
  # Plot panel a, mean viscosity at 75 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdmath $logetagrid 40 AND = junk.grd
  #gmt grdfilter $logetagrid -D1 -Fg200 -Gjunk.grd
  #logetagrid=junk.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y32c -X10c > $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 75 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=150
  # Plot panel b, mean viscosity at 150 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx-4.5c/-1.1c+w7c/0.25c+e+h -Cvisc.cpt -B1f0.25+l"log@-10@-@~\155@~@-@~\104\150@~@- (Pa s)" -Al -O -K >> $ps
  echo "0.6 1.05 150 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=250
  # Plot panel c, mean viscosity at 250 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y-12.5c -X-12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 250 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 c" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=350
  # Plot panel d, mean viscosity at 350 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "0.6 1.05 350 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 d" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps

  gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
  rm *.cpt

  # Set plot name
  f_name="${grid_type}_quad_panel_std"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  gmt makecpt -Cpolar -D -T0/1.0/0.05 -I > visc.cpt

  depth=75
  # Plot panel a, std viscosity at 75 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  logetagrid_mean=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y32c -X10c > $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  #gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Cvisc.cpt -B1f0.25+l"log@-10@-@~\150@~ (Pa s)" -Al -O -K >> $ps
  echo "-0.05 1.05 75 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=150
  # Plot panel b, std viscosity at 150 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  logetagrid_mean=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx-4.5c/-1.1c+w7c/0.25c+e+h -Cvisc.cpt -B0.25f0.25+l"log@-10@-@~\163@~@-@~\104\150@~@-" -Al -O -K >> $ps
  echo "0.6 1.05 150 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=250
  # Plot panel c, std viscosity at 250 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  logetagrid_mean=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y-12.5c -X-12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  #gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Cvisc.cpt -B1f0.25+l"log@-10@-@~\150@~ (Pa s)" -Al -O -K >> $ps
  echo "-0.05 1.05 250 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 c" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=350
  # Plot panel d, std viscosity at 350 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  logetagrid_mean=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  #gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Cvisc.cpt -B1f0.25+l"log@-10@-@~\150@~ (Pa s)" -Al -O -K >> $ps
  echo "0.6 1.05 350 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 d" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps

  gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
  rm *.cpt

}

plot_temperature_quad_panel(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  grid_type="temperature"
  full_grid_type="potential_temperature"
  fold_grid=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  fold_grid_visc=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/viscosity_grids/${summary_type}

  # Set plot name
  fold_actual_plot_output=${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  f_name="${grid_type}_quad_panel_mean"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  gmt makecpt -Cpolar -D -T1.0/1.7/0.05 > temp.cpt

  depth=75
  # Plot panel a, mean temperature at 75 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y32c -X10c > $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  logetagrid=${fold_grid_visc}/viscosity_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 75 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps

  depth=150
  # Plot panel b, mean temperature at 150 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  logetagrid=${fold_grid_visc}/viscosity_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx-4.5c/-1.1c+w7c/0.25c+e+h -Ctemp.cpt -B0.1f0.05+l"@~\155@~@-@~\161@~@- / 1000 (@~\260@~C)" -Al -O -K >> $ps
  echo "0.6 1.05 150 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps

  depth=250
  # Plot panel c, mean temperature at 250 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y-12.5c -X-12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  logetagrid=${fold_grid_visc}/viscosity_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 250 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 c" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps

  depth=350
  # Plot panel d, mean temperature at 350 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  logetagrid=${fold_grid_visc}/viscosity_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "0.6 1.05 350 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 d" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps

  gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
  rm *.cpt

  # Set plot name
  fold_actual_plot_output=${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  f_name="${grid_type}_quad_panel_std"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  gmt makecpt -Cpolar -D -T0/60/2.5 > temp.cpt

  depth=75
  # Plot panel a, std temperature at 75 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y32c -X10c > $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  logetagrid=${fold_grid_visc}/viscosity_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 75 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=150
  # Plot panel b, std temperature at 150 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  logetagrid=${fold_grid_visc}/viscosity_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx-4.5c/-1.1c+w7c/0.25c+ef+h -Ctemp.cpt -B10f5+l"@~\163@~@-@~\161@~@- (@~\260@~C)" -Al -O -K >> $ps
  echo "0.6 1.05 150 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=250
  # Plot panel c, std temperature at 250 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y-12.5c -X-12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  logetagrid=${fold_grid_visc}/viscosity_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 250 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 c" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=350
  # Plot panel d, std temperature at 350 km
  Tgrid=${fold_grid}/${full_grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $Tgrid -Ctemp.cpt $rgn_map $proj_map -O -K >> $ps
  logetagrid=${fold_grid_visc}/viscosity_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "0.6 1.05 350 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 d" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps

  gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
  rm *.cpt

}

plot_viscosity_quad_panel(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  grid_type="viscosity"
  testgridfold=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/distribution/individual
  fold_grid=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}

  # Set plot name
  fold_actual_plot_output=${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  f_name="${grid_type}_quad_panel_mean"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  gmt makecpt -Chot -D -T18/23/0.25 -I -G0.08/1 > visc.cpt

  depth=75
  # Plot panel a, mean viscosity at 75 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  #gmt grdmath $logetagrid 40 AND = junk.grd
  #gmt grdfilter $logetagrid -D1 -Fg200 -Gjunk.grd
  #logetagrid=junk.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y32c -X10c > $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 75 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=150
  # Plot panel b, mean viscosity at 150 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx-4.5c/-1.1c+w7c/0.25c+e+h -Cvisc.cpt -B1f0.25+l"log@-10@-@~\155@~@-@~\150@~@- (Pa s)" -Al -O -K >> $ps
  echo "0.6 1.05 150 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=250
  # Plot panel c, mean viscosity at 250 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y-12.5c -X-12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 1.05 250 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 c" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=350
  # Plot panel d, mean viscosity at 350 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt grdcontour $logetagrid $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "0.6 1.05 350 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 d" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps

  gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
  rm *.cpt

  # Set plot name
  fold_actual_plot_output=${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/${summary_type}
  f_name="${grid_type}_quad_panel_std"
  ps="${fold_actual_plot_output}/${f_name}.ps"
  jpg="${fold_actual_plot_output}/${f_name}.jpg"

  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  gmt makecpt -Chot -D -T0/1.0/0.05 -I -G0.08/1 > visc.cpt

  depth=75
  # Plot panel a, std viscosity at 75 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  logetagrid_mean=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y32c -X10c > $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt grdcontour $logetagrid_mean $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  #gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Cvisc.cpt -B1f0.25+l"log@-10@-@~\150@~ (Pa s)" -Al -O -K >> $ps
  echo "-0.05 1.05 75 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=150
  # Plot panel b, std viscosity at 150 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  logetagrid_mean=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt grdcontour $logetagrid_mean $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx-4.5c/-1.1c+w7c/0.25c+ef+h -Cvisc.cpt -B0.25f0.25+l"log@-10@-@~\163@~@-@~\150@~@-" -Al -O -K >> $ps
  echo "0.6 1.05 150 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=250
  # Plot panel c, std viscosity at 250 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  logetagrid_mean=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y-12.5c -X-12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt grdcontour $logetagrid_mean $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  #gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Cvisc.cpt -B1f0.25+l"log@-10@-@~\150@~ (Pa s)" -Al -O -K >> $ps
  echo "-0.05 1.05 250 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "-0.05 -0.05 c" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  depth=350
  # Plot panel d, std viscosity at 350 km
  logetagrid=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
  logetagrid_mean=${fold_grid}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K -Y0c -X12c >> $ps
  gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
  gmt grdcontour $logetagrid_mean $rgn_map $proj_map -C+22.5 -W1p,177/252/3 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  #gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Cvisc.cpt -B1f0.25+l"log@-10@-@~\150@~ (Pa s)" -Al -O -K >> $ps
  echo "0.6 1.05 350 km" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  echo "0.715 -0.05 d" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps

  gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
  rm *.cpt

}

# Load general parameters and filepaths from config file
rgn="-R0/360/-90/-60"
inc="-I0.05d/0.05d"
model="ANT-20"
fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_scripts=${fold_base}/scripts
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')
sol50=$(awk '$1 ~ /^solidus_50km/' config.ini | awk '{print $3}')
date=$(awk '$1 ~ /^date/' config.ini | awk '{print $3}')
Tp=$(awk '$1 ~ /^potential_temperature/' config.ini | awk '{print $3}')

compare_Tp=$(awk '$1 ~ /^compare_potential_temperature/' config.ini | awk '{print $3}')
compare_sol50=$(awk '$1 ~ /^compare_solidus_50km/' config.ini | awk '{print $3}')
compare_date=$(awk '$1 ~ /^compare_date/' config.ini | awk '{print $3}')

mkdir -p ${fold_plot_output}
mkdir -p ${fold_plot_output}/${date}
mkdir -p ${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}

for grid_type in density temperature viscosity; do
    mkdir -p ${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids
    mkdir -p ${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/distribution
    mkdir -p ${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/distribution/individual
    mkdir -p ${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/distribution/summary
    mkdir -p ${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}/${grid_type}_grids/MAP
done

summary_type="distribution/summary" # choose distribution/(individual/summary) or MAP

#convert_temperature_potential_temperature_grids
#plot_temperature_quad_panel
plot_viscosity_quad_panel
#plot_thermodynamic_histograms
#plot_thermodynamic_distance_histograms

#get_comparison_plot_output
#convert_temperature_potential_temperature_comparison_grids
#plot_temperature_distance_quad_panel
#plot_viscosity_distance_quad_panel