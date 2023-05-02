#!/bin/bash

set_hoggard_FR12_2020_parameter_values_YT16(){
    #specific to YT16 parameterisation
    echo 69.0 -11.7 2.83 22.117 963 0.0 4.5 > ${fold_data_output_histogram}/params_FR12_hoggard_2020.txt
}

set_hazzard_ANT-20_2022_parameter_values_YT16(){
    #specific to YT16 parameterisation
    echo 74.8 -13.1 2.09 22.9 476 5.02 1.65 > ${fold_data_output_histogram}/params_ANT-20_hazzard_2022.txt
}

set_hazzard_FR12_parameter_values_YT16(){
    #specific to YT16 parameterisation
    awk 'BEGIN {ORS=" "}; {if (NR>1 && NR<9) {print $2}}' ${fold_data_output_summary}/MAP_model.txt > junk
    echo -n > ${fold_data_output_histogram}/params_${output_time}.txt 
    for param in $(seq 1 1 7); do
        unit=$(awk '{print $'$param'}' ${fold_data_output_histogram}/param_units.txt)
        value=$(awk '{print $'$param'}' junk)
        awk 'BEGIN {ORS=" "}; {print '$unit'*$'$param'}' junk >> ${fold_data_output_histogram}/params_${output_time}.txt 
    done 
}

set_parameter_annotations_YT16(){
    #specific to YT16 parameterisation
    echo -e "@~\155@~@-0@- (GPa)\n@~\266@~@~\155@~/@~\266@~T (MPa / @~\260@~ C)\n@~\266@~@~\155@~/@~\266@~P\nlog@-10@-@~\150@~@-r@- (Pa s)\nE@-a@- (kJ / mol)\nV@-a@- (cm@+3@+ / mol)\n@~\266@~T@-s@-/@~\266@~z (@~\260@~ C / km)" > ${fold_data_output_histogram}/param_annots.txt
}

set_parameter_units_YT16(){
    #specific to YT16 parameterisation
    echo 1 1000 1 1 0.001 1000000 1 > ${fold_data_output_histogram}/param_units.txt 
}

set_parameter_ticks_YT16(){
    #specific to YT16 parameterisation
    echo 5 5 0.2 2 250 5 2 > ${fold_data_output_histogram}/param_ticks.txt
    echo 1 1 0.1 1 250 5 1 >> ${fold_data_output_histogram}/param_ticks.txt
}

set_parameter_bin_widths_YT16(){
    #specific to YT16 parameterisation
    echo 1 0.001 0.1 1 250e3 5e-6 0.2 > ${fold_data_output_histogram}/param_bin_widths.txt
}

plot_parameter_histograms(){

    # Set GMT plotting parameters   
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

    ps="${fold_plot_output_histogram}/parameter_histograms.ps"
    jpg="${fold_plot_output_histogram}/parameter_histograms.jpg"

    scale="-JX10c/10c"

    i=0
    j=$(echo $i | awk '{print $1+1}')
    x_major_tick=5
    x_minor_tick=1
    bin_width=0.25
    awk '(NR > 1) {print $'$j'}' ${file_BANCAL22_data} | gmt pshistogram -Z1 $scale -W${bin_width} -I > junk
    x_min=$(awk '{print $1}' junk)
    x_max=$(awk '{print $2}' junk)
    y_max=$(awk '{print 1.05*$4}' junk)
    rgn="-R${x_min}/${x_max}/0/${y_max}"
    x_label="log@-10@-(P)"
    y_label="f"
    gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -BWSne -X25c -Y65c -K > $ps
    awk '{print $'$i'}' ${file_fwd_models} | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gred -O -K -t50 >> $ps
    awk '(NR > 1) {print $'$j'}' ${file_BANCAL22_data} | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gblue -O -K -t50 >> $ps
    echo $i

    for i in $(seq 1 1 7); do

        echo "Plotting parameter ${i}"
        move=$(echo $i | awk '{if ($1 % 2 == 0) print "-X-11c -Y-12c"; else print "-X11c"}')
        unit=$(awk '{print $'$i'}' ${fold_data_output_histogram}/param_units.txt)
        j=$(echo $i | awk '{print $1+1}')
        mu=$(awk 'NR==2 {print $'$i'}' ${file_priors})
        sigma=$(awk 'NR==3 {print $'$i'}' ${file_priors})
        x_min=$(echo $mu $sigma $unit | awk '{print $3*($1-4*$2)}')
        x_max=$(echo $mu $sigma $unit | awk '{print $3*($1+4*$2)}')
        echo $x_min $x_max $mu $sigma
        bin_width=$(echo $x_min $x_max | awk '{print ($2-$1)/100}')
        x_major_tick=$(awk '(NR == 1) {print $'$i'}' ${fold_data_output_histogram}/param_ticks.txt)
        x_minor_tick=$(awk '(NR == 2) {print $'$i'}' ${fold_data_output_histogram}/param_ticks.txt)
        rgn="-R${x_min}/${x_max}/0/100"
        awk '(NR > 1) {print '$unit'*$'$j'}' ${file_BANCAL22_data} | gmt pshistogram -F -Z1 $rgn $scale -W${bin_width} -I > junk
        junk_x_min=$(awk '{print $1}' junk)
        junk_x_max=$(awk '{print $2}' junk)
        x_min=$(echo $x_min $junk_x_min | awk '{if ($1 < $2) print $1; else print $2}')
        x_max=$(echo $x_max $junk_x_max | awk '{if ($1 > $2) print $1; else print $2}')
        y_max=$(awk '{print 1.02*$4}' junk)
        rgn="-R${x_min}/${x_max}/0/${y_max}"
        x_label=$(awk 'NR=='$i' {print $0}' ${fold_data_output_histogram}/param_annots.txt)
        y_label="f"
        x_min_data=0
        x_max_data=2
        rgn_data="-R${x_min_data}/${x_max_data}/0/${y_max}"
        gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -BWSne $move -O -K >> $ps
        #hogg=$(awk '{print $'$i'}' ${fold_data_output_histogram}/params_FR12_hoggard_2020.txt)
        #hazz_ant=$(awk '{print $'$i'}' ${fold_data_output_histogram}/params_ANT-20_hazzard_2022.txt)
        #hazz_fr=$(awk '{print $'$i'}' ${fold_data_output_histogram}/params_FR12_hazzard_${output_time}.txt)
        #echo -e "${hogg} 0\n ${hogg} ${y_max}" > junk_hogg
        #echo -e "${hazz_ant} 0\n ${hazz_ant} ${y_max}" > junk_hazz_ant
        #echo -e "${hazz_fr} 0\n ${hazz_fr} ${y_max}" > junk_hazz_fr
        #gmt psxy $rgn $scale -W1.5p,black,-- junk_hogg -O -K >> $ps
        #gmt psxy $rgn $scale -W2.0p,black,.- junk_hazz_ant -O -K >> $ps
        #gmt psxy $rgn $scale -W2.5p,black,.. junk_hazz_fr -O -K >> $ps
        awk '{print '$unit'*$'$i'}' ${file_fwd_models} | gmt pshistogram -F -Z1 $rgn $scale -W${bin_width} -Gred -O -K -t50 >> $ps
        awk '(NR > 1) {print '$unit'*$'$j'}' ${file_BANCAL22_data} | gmt pshistogram -F -Z1 $rgn $scale -W${bin_width} -Gblue -O -K -t50 >> $ps
        echo $i
    done

    gmt psbasemap $rgn $scale -B0 -O >> $ps
    gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_BANCAL22_runs=$(awk '$1 ~ /^BANCAL22_runs/' config.ini | awk '{print $3}')
fold_fwd_model_runs=$(awk '$1 ~ /^fwd_models/' config.ini | awk '{print $3}')
fold_date=$(awk '$1 ~ /^date/' config.ini | awk '{print $3}')
Tp=$(cat ${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/potential_temperature/potential_temperature.T)
sol50=$(cat ${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/data/potential_temperature/solidus_50km_temperature.T)
file_BANCAL22_data=${fold_data_input}/${fold_BANCAL22_runs}/${fold_date}/samples_postburnin.csv
file_fwd_models=${fold_data_input}/${fold_fwd_model_runs}/800k_fwd_models.txt
file_priors=${fold_data_input}/${fold_fwd_model_runs}/priors.txt
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')
output_time=$fold_date

mkdir -p ${fold_data_output}
mkdir -p ${fold_data_output}/${output_time}
mkdir -p ${fold_data_output}/${output_time}/Tp_${Tp}_sol50_${sol50}
mkdir -p ${fold_data_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/histogram
mkdir -p ${fold_plot_output}
mkdir -p ${fold_plot_output}/${output_time}/Tp_${Tp}_sol50_${sol50}
mkdir -p ${fold_plot_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/histogram

fold_data_output_summary=${fold_data_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/parameter_summary
fold_data_output_histogram=${fold_data_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/histogram 
fold_plot_output_histogram=${fold_plot_output}/${output_time}/Tp_${Tp}_sol50_${sol50}/histogram

echo $fold_data_output_histogram
set_parameter_units_YT16
set_parameter_ticks_YT16
set_parameter_bin_widths_YT16
set_parameter_annotations_YT16
plot_parameter_histograms