#!/bin/bash

convert_parameters_YT16(){
    #specific to YT16 parameterisation
    echo 1 1000 1 1 0.001 1000000 1 > ${fold_data_output_tradeoff}/param_conversions.txt
    #convert dmudT into e-2 Gpa / K, Ea into kJ / mol, Va into cm^3 / mol
}

set_parameter_ticks_YT16(){
    #specific to YT16 parameterisation
    echo 2 4 0.2 1 250 5 0.2 > ${fold_data_output_tradeoff}/param_ticks.txt
}

set_parameter_annotations_YT16(){
    #specific to YT16 parameterisation
    echo -e "@~\155@~@-0@- (GPa)\n@~\266@~@~\155@~/@~\266@~T (MPa/@~\260@~C)\n@~\266@~@~\155@~/@~\266@~P\nlog@-10@-@~\150@~@-r@- (Pa s)\nE@-a@- (kJ/mol)\nV@-a@- (cm@+3@+/mol)\n@~\266@~T@-s@-/@~\266@~z (@~\260@~C/km)" > ${fold_data_output_tradeoff}/param_annots.txt
}

calculate_trade_off(){

    parameterlines=$(wc -l ${file_BANCAL22_data} | awk '{print $1-1}')
    no_params=7
    echo $parameterlines $no_params
    #read in (x,y,z) lines where x=parameter 1, y=parameter 2, z=1 for each parameter line

    for ((i=1; i<=$no_params; i++)); do

      for ((j=1; j<=i; j++)); do

        if [[ $i -eq $j ]]; then
          echo -n
        else
          if [ -f ${fold_data_output_tradeoff}/${j}_${i}_grid_density.txt ]; then
            echo trade_offs already calculated, using these..
          else
            sf_i=$(awk '{print $('$i')}' ${fold_data_output_tradeoff}/param_conversions.txt)
            sf_j=$(awk '{print $('$j')}' ${fold_data_output_tradeoff}/param_conversions.txt)
            echo calculating trade_off for params $i, $j
            awk 'NR>=2 {print $('$i'+1), $('$j'+1), 1}' ${file_BANCAL22_data} | awk '{print $1*'$sf_i',$2*'$sf_j',$3}' > ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt
            awk '{print $2, $1, $3}' ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt > ${fold_data_output_tradeoff}/${j}_${i}_trade_off.txt
            i_min=$(sort -nk 1 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | head -n 1 | awk '{print $1}' )
            i_max=$(sort -nk 1 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | tail -n 1 | awk '{print $1}')
            j_min=$(sort -nk 2 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | head -n 1 | awk '{print $2}')
            j_max=$(sort -nk 2 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | tail -n 1 | awk '{print $2}')
            i_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/1000}')
            j_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/1000}')

            gmt blockmean ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt -R${i_min}/${i_max}/${j_min}/${j_max} -I${i_spacing}/${j_spacing} -Sn -C > ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt
            gmt blockmean ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt -R${i_min}/${i_max}/${j_min}/${j_max} -I${i_spacing}/${j_spacing} -Sn -C | awk '{print $2, $1, $3}' > ${fold_data_output_tradeoff}/${j}_${i}_grid_density.txt
          fi
        fi

      done

    done
}

plot_trade_off(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 20p FONT_ANNOT_PRIMARY 16p PS_MEDIA a0 MAP_FRAME_TYPE plain

  # Set plot name
  ps="${fold_plot_output_tradeoff}/density_trade_off_presentation.ps"
  jpg="${fold_plot_output_tradeoff}/density_trade_off_presentation.jpg"

  no_params=7
  val=$(echo $no_params | awk '{print $0-1}')
  val4=$(echo $no_params | awk '{print $0-2}')

  cpt=prob_density.cpt
  gmt makecpt -T0.2/1.8/0.025 -Chot -D -A50 > $cpt


  for ((i=1; i<=$val; i++)); do

    val2=$(echo $i | awk '{print $0-1}')
    val3=$(echo $i | awk '{print $0+1}')

    for ((j=$val3; j<=$no_params; j++)); do


        echo plotting trade_off for params $i, $j

        i_min=$(sort -nk 1 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | head -n 1 | awk '{print $1}')
        i_max=$(sort -nk 1 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | tail -n 1 | awk '{print $1}')
        j_min=$(sort -nk 2 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | head -n 1 | awk '{print $2}')
        j_max=$(sort -nk 2 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | tail -n 1 | awk '{print $2}')
        i_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/1000}')
        j_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/1000}')

        i_axis_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/3}')
        i_axis_minor_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/6}')
        j_axis_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/5}')
        j_axis_minor_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/10}')

        i_axis_spacing=$(awk '{print $('$i')}' ${fold_data_output_tradeoff}/param_ticks.txt)
        i_axis_minor_spacing=$(echo $i_axis_spacing | awk '{print $0/2}')
        j_axis_spacing=$(awk '{print $('$j')}' ${fold_data_output_tradeoff}/param_ticks.txt)
        j_axis_minor_spacing=$(echo $j_axis_spacing | awk '{print $0/2}')

        k_min=$(sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | head -n 1 | awk '{print log($3)/log(10)}')
        k_max=$(sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | tail -n 1 | awk '{print log($3)/log(10)}')
        k_spacing=$(echo $k_min $k_max | awk '{print -($1-$2)/100}')

        rgn="-R${i_min}/${i_max}/${j_min}/${j_max}"
        scale="-JX5.3c/5.3c"
        plot_dist=$(echo 5.5 | awk '{print $1}')

        xval=$(echo $i $plot_dist | awk '{print $1*$2+1}')
        yval=$(echo $j $plot_dist | awk '{print 80-$1*$2}')

        i_annot=$(awk 'NR=='$i' {print $0}' ${fold_data_output_tradeoff}/param_annots.txt)
        j_annot=$(awk 'NR=='$j' {print $0}' ${fold_data_output_tradeoff}/param_annots.txt)


        if [[ $i -eq 1 ]] && [[ $j -eq 2 ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BWrbt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -K > $ps
        elif [[ $i -eq 1 ]] && [[ $j -ge 3 ]] && [[ $j -le $val ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BWrbt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        elif [[ $i -eq 1 ]] && [[ $j -eq $no_params ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BWrSt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        elif [[ $j -eq $no_params ]] && [[ $i -ge 2 ]] && [[ $i -le $val4 ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BlrSt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        elif [[ $i -eq $val ]] && [[ $j -eq $no_params ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BlrSt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        elif [[ $i -ge 2 ]] && [[ $i -le $val4 ]] && [[ $j -eq $val3 ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -Blrbt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        else
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -Blrbt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        fi

        sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | awk '{print $1, $2, log($3)/log(10)}' | gmt psxy -Xa${xval}c -Ya${yval}c $rgn $scale -Sc0.1c -C$cpt -O -K >> $ps

        if [[ $i -eq $val ]] && [[ $j -eq $no_params ]]; then
          gmt psscale -C$cpt -Dx2c/51c+w15c/0.55c+e+m -R -J -O -Bxa0.5+l"log@-10@-(@~\162@~@-sample@-)" >> $ps
        fi
        rm -f junk

    done
  done
	gmt psconvert $ps -Tj -A0.25c -P -Z -E600
  rm -f gmt*

}

plot_trade_off_standard_range(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 20p FONT_ANNOT_PRIMARY 16p PS_MEDIA a0 MAP_FRAME_TYPE plain

  # Set plot name
  ps="${fold_plot_output_tradeoff}/density_trade_off_standard_range.ps"
  jpg="${fold_plot_output_tradeoff}/density_trade_off_standard_range.jpg"

  no_params=7
  val=$(echo $no_params | awk '{print $0-1}')
  val4=$(echo $no_params | awk '{print $0-2}')

  cpt=prob_density.cpt
  gmt makecpt -T0.2/1.8/0.025 -Chot -D -A50 > $cpt


  for ((i=1; i<=$val; i++)); do

    val2=$(echo $i | awk '{print $0-1}')
    val3=$(echo $i | awk '{print $0+1}')

    for ((j=$val3; j<=$no_params; j++)); do


        echo plotting trade_off for params $i, $j

        n_outputs=$(ls -d ${fold_data_output}/*/ | wc -l)

        for ((k=1; k<=$n_outputs; k++)); do

          tmp_fold=$(ls -ld ${fold_data_output}/*/ | awk '{if (NR=='${k}') print $9}')
          file=${tmp_fold}${i}_${j}_trade_off.txt

          if [ $k -eq 1 ]; then

            i_min=$(sort -nk 1 ${file} | head -n 1 | awk '{print $1}')
            i_max=$(sort -nk 1 ${file} | tail -n 1 | awk '{print $1}')
            j_min=$(sort -nk 2 ${file} | head -n 1 | awk '{print $2}')
            j_max=$(sort -nk 2 ${file} | tail -n 1 | awk '{print $2}')

          else

            tmp_i_min=$(sort -nk 1 ${file} | head -n 1 | awk '{print $1}')
            tmp_i_max=$(sort -nk 1 ${file} | tail -n 1 | awk '{print $1}')
            tmp_j_min=$(sort -nk 2 ${file} | head -n 1 | awk '{print $2}')
            tmp_j_max=$(sort -nk 2 ${file} | tail -n 1 | awk '{print $2}')
            i_min=$(echo $i_min $tmp_i_min | awk '{if ($2 < $1) print $2; else print $1}')
            i_max=$(echo $i_max $tmp_i_max | awk '{if ($2 > $1) print $2; else print $1}')
            j_min=$(echo $j_min $tmp_j_min | awk '{if ($2 < $1) print $2; else print $1}')
            j_max=$(echo $j_max $tmp_j_max | awk '{if ($2 > $1) print $2; else print $1}')

          fi

        done

        i_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/1000}')
        j_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/1000}')

        i_axis_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/3}')
        i_axis_minor_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/6}')
        j_axis_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/5}')
        j_axis_minor_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/10}')

        i_axis_spacing=$(awk '{print $('$i')}' ${fold_data_output_tradeoff}/param_ticks.txt)
        i_axis_minor_spacing=$(echo $i_axis_spacing | awk '{print $0/2}')
        j_axis_spacing=$(awk '{print $('$j')}' ${fold_data_output_tradeoff}/param_ticks.txt)
        j_axis_minor_spacing=$(echo $j_axis_spacing | awk '{print $0/2}')

        k_min=$(sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | head -n 1 | awk '{print log($3)/log(10)}')
        k_max=$(sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | tail -n 1 | awk '{print log($3)/log(10)}')
        k_spacing=$(echo $k_min $k_max | awk '{print -($1-$2)/100}')

        rgn="-R${i_min}/${i_max}/${j_min}/${j_max}"
        scale="-JX5.3c/5.3c"
        plot_dist=$(echo 5.5 | awk '{print $1}')

        xval=$(echo $i $plot_dist | awk '{print $1*$2+1}')
        yval=$(echo $j $plot_dist | awk '{print 80-$1*$2}')

        i_annot=$(awk 'NR=='$i' {print $0}' ${fold_data_output_tradeoff}/param_annots.txt)
        j_annot=$(awk 'NR=='$j' {print $0}' ${fold_data_output_tradeoff}/param_annots.txt)


        if [[ $i -eq 1 ]] && [[ $j -eq 2 ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BWrbt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -K > $ps
        elif [[ $i -eq 1 ]] && [[ $j -ge 3 ]] && [[ $j -le $val ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BWrbt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        elif [[ $i -eq 1 ]] && [[ $j -eq $no_params ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BWrSt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        elif [[ $j -eq $no_params ]] && [[ $i -ge 2 ]] && [[ $i -le $val4 ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BlrSt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        elif [[ $i -eq $val ]] && [[ $j -eq $no_params ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BlrSt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        elif [[ $i -ge 2 ]] && [[ $i -le $val4 ]] && [[ $j -eq $val3 ]]; then
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -Blrbt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        else
          gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -Blrbt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -O -K >> $ps
        fi

        sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | awk '{print $1, $2, log($3)/log(10)}' | gmt psxy -Xa${xval}c -Ya${yval}c $rgn $scale -Sc0.1c -C$cpt -O -K >> $ps

        if [[ $i -eq $val ]] && [[ $j -eq $no_params ]]; then
          gmt psscale -C$cpt -Dx2c/51c+w15c/0.55c+e+m -R -J -O -Bxa0.5+l"log@-10@-(@~\162@~@-sample@-)" >> $ps
        fi
        rm -f junk

    done
  done
	gmt psconvert $ps -Tj -A0.25c -P -Z -E600
  rm -f gmt*

}

plot_trade_off_single(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 12p FONT_ANNOT_PRIMARY 12p PS_MEDIA a0 MAP_FRAME_TYPE plain

  i=4
  j=5

  # Set plot name
  ps="${fold_plot_output_tradeoff}/density_trade_off_presentation_${i}_${j}.ps"
  jpg="${fold_plot_output_tradeoff}/density_trade_off_presentation_${i}_${j}.jpg"

  no_params=7
  val=$(echo $no_params | awk '{print $0-1}')
  val4=$(echo $no_params | awk '{print $0-2}')

  cpt=prob_density.cpt
  gmt makecpt -T0.2/1.8/0.025 -Chot -D -A50 > $cpt

        echo plotting trade_off for params $i, $j

        i_min=$(sort -nk 1 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | head -n 1 | awk '{print $1}')
        i_max=$(sort -nk 1 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | tail -n 1 | awk '{print $1}')
        j_min=$(sort -nk 2 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | head -n 1 | awk '{print $2}')
        j_max=$(sort -nk 2 ${fold_data_output_tradeoff}/${i}_${j}_trade_off.txt | tail -n 1 | awk '{print $2}')
        i_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/1000}')
        j_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/1000}')

        i_axis_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/3}')
        i_axis_minor_spacing=$(echo $i_min $i_max | awk '{print -($1-$2)/6}')
        j_axis_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/5}')
        j_axis_minor_spacing=$(echo $j_min $j_max | awk '{print -($1-$2)/10}')

        i_axis_spacing=$(awk '{print $('$i')}' ${fold_data_output_tradeoff}/param_ticks.txt)
        i_axis_minor_spacing=$(echo $i_axis_spacing | awk '{print $0/2}')
        j_axis_spacing=$(awk '{print $('$j')}' ${fold_data_output_tradeoff}/param_ticks.txt)
        j_axis_minor_spacing=$(echo $j_axis_spacing | awk '{print $0/2}')

        k_min=$(sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | head -n 1 | awk '{print log($3)/log(10)}')
        k_max=$(sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | tail -n 1 | awk '{print log($3)/log(10)}')
        #k_min=$(sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | head -n 1 | awk '{print $3}')
        #k_max=$(sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | tail -n 1 | awk '{print $3}')
        k_spacing=$(echo $k_min $k_max | awk '{print -($1-$2)/100}')
        #echo $i $j $k_min $k_max
        #cpt=prob_density.cpt
        #gmt makecpt -T${k_min}/${k_max}/${k_spacing} -Chot -D -A50 > $cpt

        rgn="-R${i_min}/${i_max}/${j_min}/${j_max}"
        scale="-JX5.3c/5.3c"
        plot_dist=$(echo 5.5 | awk '{print $1}')

        xval=$(echo $i $plot_dist | awk '{print $1*$2+1}')
        yval=$(echo $j $plot_dist | awk '{print 80-$1*$2}')

        i_annot=$(awk 'NR=='$i' {print $0}' ${fold_data_output_tradeoff}/param_annots.txt)
        j_annot=$(awk 'NR=='$j' {print $0}' ${fold_data_output_tradeoff}/param_annots.txt)
        #echo "$i_annot" "$j_annot"

        gmt psbasemap -Xa${xval}c -Ya${yval}c $rgn $scale -Ba0 -BWrSt -Bx${i_axis_spacing}f${i_axis_minor_spacing}+l"$i_annot" -By${j_axis_spacing}f${j_axis_minor_spacing}+l"$j_annot" -K > $ps
        sort -nk 3 ${fold_data_output_tradeoff}/${i}_${j}_grid_density.txt | awk '{print $1, $2, log($3)/log(10)}' | gmt psxy -Xa${xval}c -Ya${yval}c $rgn $scale -Sc0.07c -C$cpt -O >> $ps
        rm -f junk

	gmt psconvert $ps -Tj -A0.05c -P -Z -E600
  rm -f gmt*

}

fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_BANCAL22_runs=$(awk '$1 ~ /^BANCAL22_runs/' config.ini | awk '{print $3}')
fold_date=$(awk '$1 ~ /^date/' config.ini | awk '{print $3}')
Tp=$(cat ${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/potential_temperature/potential_temperature.T)
sol50=$(cat ${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/potential_temperature/solidus_50km_temperature.T)
file_BANCAL22_data=${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/samples_postburnin.csv
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')
output_time=$fold_date

mkdir -p ${fold_data_output}
mkdir -p ${fold_data_output}/${output_time}
mkdir -p ${fold_data_output}/${output_time}/Tp_${Tp}_sol50_${sol50}
mkdir -p ${fold_data_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/trade_off
mkdir -p ${fold_plot_output}
mkdir -p ${fold_plot_output}/${output_time}/Tp_${Tp}_sol50_${sol50}
mkdir -p ${fold_plot_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/trade_off

fold_data_output_tradeoff=${fold_data_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/trade_off
fold_plot_output_tradeoff=${fold_plot_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/trade_off

convert_parameters_YT16
set_parameter_ticks_YT16
set_parameter_annotations_YT16
calculate_trade_off
plot_trade_off
#plot_trade_off_standard_range

rm -f junk* gmt.* *.cpt