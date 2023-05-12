#!/bin/bash

plot_sample_locs(){

    # Set GMT plotting parameters
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

    plot_name="${fold_plot_output}/sample_locs"
    ps="${plot_name}.ps"
    jpg="${plot_name}.jpg"

    proj_map="-JS0/-90/10c"
    rgn_map="-R0/360/-90/-60"
    rgnx="-R0/1/0/1"
    scalex="-JX15c/10c"

    f_sample_locs=${fold_data_input}/${input_locs}/sample_locs.txt
    no_locs=$(wc -l ${f_sample_locs} | awk '{print $1}')

    gmt psbasemap $rgn_map $proj_map -B0 -X20c -Y20c -K > $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps

    # loop over each location we wish to sample, extracting coordinates
    for ((i=1; i<=${no_locs}; i++)); do

        j=$(echo $i | awk '{print $1+1}')
        loc=$(awk '{if (NR=='${j}') print $1}' ${f_sample_locs})
        lon=$(awk '{if (NR=='${j}') print $2}' ${f_sample_locs})
        lat=$(awk '{if (NR=='${j}') print $3}' ${f_sample_locs})

        echo ${lon} ${lat} 100k | gmt psxy $rgn_map $proj_map -SE- -Gblue -O -K >> $ps

    done

    gmt psbasemap $rgn_map $proj_map -Bxa30f15g5 -Bya10f5g5 -O >> $ps
    gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
    
}

fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
input_locs=$(awk '$1 ~ /^input_locations/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')

fold_plot_output=${fold_plot_output}/sample_locs
mkdir -p ${fold_plot_output}

plot_sample_locs