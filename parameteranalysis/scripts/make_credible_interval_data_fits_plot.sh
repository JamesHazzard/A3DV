#!/bin/bash
#source ~/.bash_profile

plot_credible_interval(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0

  mkdir -p ${fold_plot_output_cidf}

  # Set plot name
  ps="${fold_plot_output_cidf}/CI_data_fit.ps"
  jpg="${fold_plot_output_cidf}/CI_data_fit.jpg"

  # Set colour array for plot
  col_arr=(black black black)
  envelope_col_arr=(blue purple palevioletred)

  ##############################################################################################################
  # Vs-T plot
  ##############################################################################################################
  rgn="-R775/1475/4.15/4.7"
  scale="-JX21.5c/13.4375c"
  gmt psbasemap $rgn $scale -B0 -K -Y25c > $ps
  depthlines=$(wc -l depth_slices | awk '{print $1}')

  i=0
  j=${i}
  for ((l=1; l<=$depthlines; l++)); do

    z=$(awk 'NR=='$l'{printf "%.1f", $1}' depth_slices)
    z_min=$(echo $z | awk '{printf "%.0f", $1-12.5}')
    z_max=$(echo $z | awk '{printf "%.0f", $1+12.5}')
    #awk '{if ($2=='$z'){print $4, $3, $5}}' ${inputfold}/plate_${model}${end}.txt | \
    awk '{if ($4=='$z'){print $3, $1, $2}}' ${fold_BANCAL22_input_data}/plate/plate.VseTz | \
        gmt psxy $rgn $scale -G${col_arr[$i]} -Ey1p/0.5p -K -O >> $ps
    awk '{print $1,$2,$3,$4}' ${fold_data_output_cidf_envelope}/CI_99_envelope_plate_${z_min}_${z_max}.TVs |  \
    gmt psxy $rgn $scale -G${envelope_col_arr[$i]} -O -K -L+b+t10 -t50 >> $ps
    awk '{print $1,$2,$3,$4}' ${fold_data_output_cidf_envelope}/CI_50_envelope_plate_${z_min}_${z_max}.TVs | \
    gmt psxy $rgn $scale -G${envelope_col_arr[$j]} -O -K -L+b+t10 -t40 >> $ps
    #awk '{if ($2=='$z'){print $4, $3, $5}}' ${inputfold}/plate_${model}${end}.txt | \
    awk '{if ($4=='$z'){print $3, $1, $2}}' ${fold_BANCAL22_input_data}/plate/plate.VseTz | \
        gmt psxy $rgn $scale -G${col_arr[$i]} -Sc0.225c -Ey0p/0p -K -O >> $ps
    i=${i}+1
    j=${i}

  done

  echo "780.5 4.692 a" | gmt pstext $rgn $scale -F+f18p+jTL -Gwhite -TO -W1p -O -K >> $ps
  gmt psbasemap $rgn $scale -Bpx100f50+l"Temperature (@~\260@~C)" -Bpy0.2f0.1+l"V@-S@- (km s@+-1@+)" -BWsNe -O -K >> $ps

  echo Finished plate model plotting

  ###############################################################################################################
  ## Adiabat plot
  ###############################################################################################################
  mindep=225
  rgn="-R4.3/5.0/145/405"
  rgnx="-R0/1/0/1"
  scale="-JX6.5c/-6.5c"
  scalex="-JX6.5c/6.5c"
  gmt psbasemap $rgn $scale -B0 -O -K -Y-7.5c >> $ps

  awk '{print $2, $4}' ${fold_data_output_cidf_envelope}/CI_99_envelope_adiabat.Vsz > temp_adiabat
  awk '{print $3, $4}' ${fold_data_output_cidf_envelope}/CI_99_envelope_adiabat.Vsz | tac >> temp_adiabat
  head -n 1 temp_adiabat >> temp_adiabat
  awk '{print $1, $2}' temp_adiabat | gmt psxy $rgn $scale -Gpalegreen2 -O -K -t50 >> $ps
  rm temp_adiabat
  awk '{print $2, $4}' ${fold_data_output_cidf_envelope}/CI_50_envelope_adiabat.Vsz > temp_adiabat
  awk '{print $3, $4}' ${fold_data_output_cidf_envelope}/CI_50_envelope_adiabat.Vsz | tac >> temp_adiabat
  head -n 1 temp_adiabat >> temp_adiabat
  awk '{print $1, $2}' temp_adiabat | gmt psxy $rgn $scale -Gpalegreen2 -O -K -t40 >> $ps
  #awk '{print $1, $2, $4}' ${inputfold}/adiabat_${model}.txt | gmt psxy $rgn $scale -Sc0.225c \
  awk '{print $1, $4, $2}' ${fold_BANCAL22_input_data}/adiabat/adiabat.VseTz | gmt psxy $rgn $scale -Sc0.225c \
     -Gblack -K -O -N -By50f25+l"Depth (km)" -Bx0.2f0.1+l"V@-S@- (km s@+-1@+)" -BWSne -Ex2p/1p >> $ps

  echo "0.025 0.97 b" | gmt pstext $rgnx $scalex -F+f18p+jTL -Gwhite -W1p -TO -O -K >> $ps
  rm -f envelope out*.temp temp_adiabat

  echo Finished adiabatic model plotting

  ###############################################################################################################
  ## Attenuation Plot
  ###############################################################################################################
  rgn="-R1e-3/1e-1/145/405"
  scale="-JX6.5cl/-6.5c"
  gmt psbasemap $rgn $scale -B0 -K -O -X7.5c >> $ps


  awk '{print $2, $4}' ${fold_data_output_cidf_envelope}/CI_99_envelope_attenuation.Qz > temp_attenuation
  awk '{print $3, $4}' ${fold_data_output_cidf_envelope}/CI_99_envelope_attenuation.Qz | tac >> temp_attenuation
  head -n 1 temp_attenuation >> temp_attenuation
  awk '{print $1, $2}' temp_attenuation | gmt psxy $rgn $scale -Gpalegreen2 -O -K -t50 >> $ps
  rm temp_attenuation
  awk '{print $2, $4}' ${fold_data_output_cidf_envelope}/CI_50_envelope_attenuation.Qz > temp_attenuation
  awk '{print $3, $4}' ${fold_data_output_cidf_envelope}/CI_50_envelope_attenuation.Qz | tac >> temp_attenuation
  head -n 1 temp_attenuation >> temp_attenuation
  awk '{print $1, $2}' temp_attenuation | gmt psxy $rgn $scale -Gpalegreen2 -O -K -t40 >> $ps
  #awk '{if ($3>=150 && $3%25==0 && $3<=400){print $1, $3, $5}}' ${inputfold}/attenuation_${model}.txt \
  awk '{if ($4>=150 && $4%25==0 && $4<=400){print $1, $4, $2}}' ${fold_BANCAL22_input_data}/attenuation/attenuation.QeVsz \
     | gmt psxy $rgn $scale -Sc0.225c -Gblack -K -O -Bx1f3p+l"Attenuation (Q@-S@-@+-1@+)" -By50f25+l"Depth (km)" -BwSne -N -Ex2p/1p >> $ps

  echo "0.025 0.97 c" | gmt pstext $rgnx $scalex -F+f18p+jTL -Gwhite -W1p -TO -O -K >> $ps
  rm -f envelope out*.temp

  echo Finished attenuation model plotting

  ###############################################################################################################
  ## Viscosity Plot
  ###############################################################################################################
  rgn="-R18/22/0/1"
  scale="-JX6.5c/-6.5c"
  gmt psbasemap $rgn $scale -B0 -K -O -X7.5c >> $ps

  v_min=$(awk 'NR==2 {print $1}' ${fold_data_output_cidf_envelope}/CI_99_envelope_viscosity.vz)
  v_max=$(awk 'NR==3 {print $1}' ${fold_data_output_cidf_envelope}/CI_99_envelope_viscosity.vz)
  z_min=0.4
  z_max=0.6
  echo $v_min $z_max > temp_viscosity
  echo $v_min $z_min >> temp_viscosity
  echo $v_max $z_min >> temp_viscosity
  echo $v_max $z_max >> temp_viscosity
  head -n 1 temp_viscosity >> temp_viscosity
  awk '{print $1, $2}' temp_viscosity | gmt psxy $rgn $scale -Gpalegreen2 -O -K -t50 >> $ps
  rm temp_viscosity
  v_min=$(awk 'NR==2 {print $1}' ${fold_data_output_cidf_envelope}/CI_50_envelope_viscosity.vz)
  v_max=$(awk 'NR==3 {print $1}' ${fold_data_output_cidf_envelope}/CI_50_envelope_viscosity.vz)
  z_min=0.4
  z_max=0.6
  echo $v_min $z_max > temp_viscosity
  echo $v_min $z_min >> temp_viscosity
  echo $v_max $z_min >> temp_viscosity
  echo $v_max $z_max >> temp_viscosity
  head -n 1 temp_viscosity >> temp_viscosity
  awk '{print $1, $2}' temp_viscosity | gmt psxy $rgn $scale -Gpalegreen2 -O -K -t40 >> $ps
  #awk '{print 20, 0.5, 1}' ${inputfold}/viscosity_${model}.txt | gmt psxy $rgn $scale -SC0.225c -Gblack -Ex2p/1p -O -K >> $ps
  awk '{print 20, 0.5, 1}' ${fold_BANCAL22_input_data}/viscosity/viscosity.neVsz | gmt psxy $rgn $scale -SC0.225c -Gblack -Ex2p/1p -O -K >> $ps

  gmt psbasemap $rgn $scale -Bpx1f0.1+l"log@-10@-@~\150@~@-UM@- (Pa s)" -BwSne -O -K >> $ps
  echo "0.025 0.97 d" | gmt pstext $rgnx $scalex -F+f18p+jTL -Gwhite -W1p -TO -O >> $ps

  echo Finished viscosity model plotting

  ################################################################################################################
  #gmt psscale -Dx+8c/+8.25c+w12c/0.25c+e -Cposterior.cpt -B50f25+l"Posterior probability" -Al -O >> $ps	#move scale to cover all 4 plots
  gmt psconvert $ps -Tj -E600 -A0.25c -P -Z

  # Clean up
  rm -f *_line* line_* *.cpt depth_slices big_points* lonlat box* intfile junk* envelope out*.temp temp_plate temp_adibat temp_attenuation temp_viscosity
}

# General parameters and filepaths
model=ANT-20
fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_scripts=$(awk '$1 ~ /^scripts/' config.ini | awk '{print $3}')
fold_conversion=${fold_base}/${fold_scripts}/vs_to_thermo_conversions
fold_scripts=${fold_base}/${fold_scripts}
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_BANCAL22_runs=$(awk '$1 ~ /^BANCAL22_runs/' config.ini | awk '{print $3}')
fold_fwd_model_runs=$(awk '$1 ~ /^fwd_models/' config.ini | awk '{print $3}')
fold_date=$(awk '$1 ~ /^date/' config.ini | awk '{print $3}')
Tp=$(cat ${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/potential_temperature/potential_temperature.T)
sol50=$(cat ${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/potential_temperature/solidus_50km_temperature.T)
file_plate=${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/plate/plate.VseTz
file_adiabat=${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/adiabat/adiabat.VseTz
file_attenuation=${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/attenuation/attenuation.QeVsz
file_viscosity=${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/viscosity/viscosity.neVsz
fold_BANCAL22_input_data=${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data
file_BANCAL22_data=${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/samples_postburnin.csv
file_fwd_models=${fold_data_input}/${fold_fwd_model_runs}/800k_fwd_models.txt
file_priors=${fold_data_input}/${fold_fwd_model_runs}/priors.txt
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')
output_time=$fold_date

fold_data_output_cidf=${fold_data_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/credible_interval_data_fits
fold_data_output_cidf_envelope=${fold_data_output_cidf}/envelope
fold_plot_output_cidf=${fold_plot_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/credible_interval_data_fits

# Mid-point depths of Vs and T slices
awk 'NR>1{a[$4]++} END{for(b in a) print b}' $file_plate | sort > depth_slices

# Top and bottom depth for Vs vs. T files
ztop=$(head -n 1 depth_slices | awk '{print int($1 - 12.5)}')
zbase=$(tail -n 1 depth_slices | awk '{print int($1 + 12.5)}')
end="_avg_noTavg_${ztop}_${zbase}"

# Run program
plot_credible_interval

rm -f junk* gmt.*