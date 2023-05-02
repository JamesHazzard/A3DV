#!/bin/bash

calculate_LAB_comparison(){

    f_xyz_pappa=${fold_data_input_LAB}/pappa19/model3-LAB_lonadjust.xyz
    f_grid_pappa=${fold_data_input_LAB}/pappa19/model3-LAB.grd
    f_grid_ANT=${fold_data_input_LAB}/ANT-20/LAB1200_mean_filtered.grd
    f_grid_SL2013sv=${fold_data_input_LAB}/SL2013sv/mean.grd
    f_grid_SL2013sv_ant=${fold_data_input_LAB}/SL2013sv/ant_mean.grd
    f_diff_ANT_pappa=${fold_data_input_LAB}/difference/ANT-20_minus_pappa19.grd
    f_diff_ANT_SL2013sv=${fold_data_input_LAB}/difference/ANT-20_minus_SL2013sv.grd

    if [[ ! -e ${f_grid_pappa} ]]; then

        gmt sphinterpolate ${f_xyz_pappa} -R0/360/-90/-60 -I30m -Q1 -T0.3 -G${f_grid_pappa}

    fi

    if [[ ! -e ${f_grid_SL2013sv_ant} ]]; then

        gmt grd2xyz ${f_grid_SL2013sv} > junk
        gmt xyz2grd -R0/360/-90/-60 -I0.5d/0.5d junk -G${f_grid_SL2013sv_ant}

    fi

    if [[ ! -e ${f_diff_ANT_pappa} ]]; then

        gmt grdmath ${f_grid_ANT} ${f_grid_pappa} SUB = ${f_diff_ANT_pappa}

    fi

    if [[ ! -e ${f_diff_ANT_SL2013sv} ]]; then

        gmt grdmath ${f_grid_ANT} ${f_grid_SL2013sv_ant} SUB = ${f_diff_ANT_SL2013sv}

    fi

    file_preamble=${f_grid_ANT}
    file_preamble_out=${fold_data_input_LAB}/ANT-20/LAB1200_mean_filtered
    gmt grd2xyz ${file_preamble} | gmt select -F${fold_data_input_polygons}/ant_polygonfile.xy > ${file_preamble_out}_ant.xyz
    python3 summary_calculator.py -a ${file_preamble_out}_ant.xyz

    file_preamble=${f_grid_pappa}
    file_preamble_out=${fold_data_input_LAB}/pappa19/model3-LAB
    gmt grd2xyz ${file_preamble} | gmt select -F${fold_data_input_polygons}/ant_polygonfile.xy > ${file_preamble_out}_ant.xyz
    python3 summary_calculator.py -a ${file_preamble_out}_ant.xyz

    file_preamble=${f_grid_SL2013sv}
    file_preamble_out=${fold_data_input_LAB}/SL2013sv/ant_mean
    gmt grd2xyz ${file_preamble} | gmt select -F${fold_data_input_polygons}/ant_polygonfile.xy > ${file_preamble_out}_ant.xyz
    python3 summary_calculator.py -a ${file_preamble_out}_ant.xyz

    file_preamble=${f_diff_ANT_pappa}
    file_preamble_out=${fold_data_input_LAB}/difference/ANT-20_minus_pappa19
    gmt grd2xyz ${file_preamble} | gmt select -F${fold_data_input_polygons}/west_ant_polygonfile.xy > ${file_preamble_out}_west.xyz
    gmt grd2xyz ${file_preamble} | gmt select -F${fold_data_input_polygons}/east_ant_polygonfile.xy > ${file_preamble_out}_east.xyz
    python3 summary_calculator.py -a ${file_preamble_out}_west.xyz
    python3 summary_calculator.py -a ${file_preamble_out}_east.xyz

    file_preamble=${f_diff_ANT_SL2013sv}
    file_preamble_out=${fold_data_input_LAB}/difference/ANT-20_minus_SL2013sv
    gmt grd2xyz ${file_preamble} | gmt select -F${fold_data_input_polygons}/west_ant_polygonfile.xy > ${file_preamble_out}_west.xyz
    gmt grd2xyz ${file_preamble} | gmt select -F${fold_data_input_polygons}/east_ant_polygonfile.xy > ${file_preamble_out}_east.xyz
    python3 summary_calculator.py -a ${file_preamble_out}_west.xyz
    python3 summary_calculator.py -a ${file_preamble_out}_east.xyz

}

plot_LAB_comparison(){

    gmt makecpt -T0/20/1 -Ccopper -I > age.cpt

    # Set GMT plotting parameters
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

    # Set plot name
    ps="${fold_plot_output}/LAB_comparison.ps"
    jpg="${fold_plot_output}/LAB_comparison.jpg"

    # Set plot area
    proj_map="-JS0/-90/10c"
    rgn_map="-R0/360/-90/-60"
    rgnx="-R0/1/0/1"
    scalex="-JX15c/10c"
    proj=$proj_map
    rgn=$rgn_map

    # Plot panel a, mean LAB ANT-20
    model="ANT-20"
    LABgrid="${fold_data_input_LAB}/ANT-20/LAB1200_mean_filtered.grd"
    gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y30c -X10c > $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    gmt psscale -Dx-1.3c/0.5c+w9c/0.25c+e+mal -Clitho.cpt -B50f10+l"z@-LAB@- (km)" -O -K >> $ps
    echo "-0.06 1.05 a" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    LABgrid_mean=$LABgrid
    LABgrid="${fold_data_input_LAB}/SL2013sv/mean.grd"
    gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y0.0c -X11.8c >> $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    echo "-0.06 1.05 b" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    model="SL2013"
    LABgrid="${fold_data_input_LAB}/ANT-20/LAB1200_std_filtered.grd"
    gmt makecpt -T-100/100/10 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y-12.3c -X-11.8c >> $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    gmt psscale -Dx-1.3c/0.5c+w9c/0.25c+e+mal -Clitho.cpt -B20f10+l"@~\104@~@-LAB@- (km)" -O -K >> $ps
    echo "-0.06 1.05 d" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    model="(ANT-20, SL2013sv)"
    LABgrid="${fold_data_input_LAB}/difference/ANT-20_minus_SL2013sv.grd"
    gmt makecpt -T-100/100/10 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y0.0c -X11.8c >> $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    echo "-0.06 1.05 e" | gmt pstext $rgnx $scalex -F+f18p+jBR -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    model="SL2013"
    LABgrid="${fold_data_input_LAB}/pappa19/model3-LAB.grd"
    gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y12.3c -X11.8c >> $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    echo "-0.06 1.05 c" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    model="(ANT-20, SL2013sv)"
    LABgrid="${fold_data_input_LAB}/difference/ANT-20_minus_pappa19.grd"
    gmt makecpt -T-100/100/10 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y-12.3c -X0c >> $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    echo "-0.06 1.05 f" | gmt pstext $rgnx $scalex -F+f18p+jBR -C+tO -Gwhite -W1p,black -N -O >> $ps
    rm litho.cpt

    gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
}

fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
input_LAB=$(awk '$1 ~ /^input_LAB/' config.ini | awk '{print $3}')
fold_data_input_LAB=${fold_data_input}/${input_LAB}
fold_data_input_polygons=${fold_data_input}/polygons    # get polygons from westvseast runs
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}') # create log file
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')

calculate_LAB_comparison
#plot_LAB_comparison