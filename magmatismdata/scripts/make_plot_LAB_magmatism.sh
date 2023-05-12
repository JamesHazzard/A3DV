#!/bin/bash

plot_LAB(){

    gmt makecpt -T0/20/1 -Ccopper -I > age.cpt

    # Set GMT plotting parameters
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

    # Set plot name
    ps="${fold_plot_output}/LAB.ps"
    jpg="${fold_plot_output}/LAB.jpg"

    # Set plot area
    proj_map="-JS0/-90/10c"
    rgn_map="-R0/360/-90/-60"
    rgnx="-R0/1/0/1"
    scalex="-JX15c/10c"
    proj=$proj_map
    rgn=$rgn_map

    # Plot panel a, mean LAB ANT-20
    model="ANT-20"
    LABgrid="${fold_input_LAB}/ANT-20/mean.grd"
    gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y30c -X10c > $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    gmt psscale -Dx0.5c/-1.8c+w9c/0.25c+ef+h -Clitho.cpt -B50f10+l"z@-LAB@- (km)" -O -K >> $ps
    echo "0.03 1.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    # Plot panel b, std LAB ANT-20
    model="SL2013"
    LABgrid="${fold_input_LAB}/ANT-20/std.grd"
    gmt makecpt -T0/50/2 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -X11.8c >> $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    gmt psscale -Dx0.5c/-1.8c+w9c/0.25c+ef+h -Clitho.cpt -B10f5+l"@~\163@~@-LAB@- (km)" -O -K >> $ps
    echo "0.63 1.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps
    rm litho.cpt

    gmt psconvert -Tj -E600 -A0.1c -P -Z $ps

}

plot_LAB_magmatism_paper_figure(){

    gmt makecpt -T0/20/1 -Ccopper -I > age.cpt

    # Set GMT plotting parameters
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

    # Set plot name
    ps="${fold_plot_output}/LAB_paper.ps"
    jpg="${fold_plot_output}/LAB_paper.jpg"

    # Set plot area
    proj_map="-JS0/-90/10c"
    rgn_map="-R0/360/-90/-60"
    rgnx="-R0/1/0/1"
    scalex="-JX15c/10c"
    proj=$proj_map
    rgn=$rgn_map

    # Plot panel a, mean LAB ANT-20
    model="ANT-20"
    LABgrid="${fold_input_LAB}/ANT-20/mean.grd"
    gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -Y30c -X10c > $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    #gmt grdcontour $LABgrid $rgn_map $proj_map -C25 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    awk '{if ($4<=10) print $1, $2, $3}' ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean.xyaszs |\
      gmt psxy $proj $rgn -Cage.cpt -St0.4c -O -K >> $ps
    gmt psscale -Dx-1.3c/0.5c+w9c/0.25c+ef+mal -Clitho.cpt -B50f10+l"z@-LAB@- (km)" -O -K >> $ps
    echo "-0.08 1.05 ANT-20" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    echo "-0.06 -0.05 a" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    # Plot panel b, std LAB
    LABgrid_mean=$LABgrid
    LABgrid="${fold_input_LAB}/SL2013/mean.grd"
    gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
    #gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y-12.3c -X0c >> $ps
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y0.0c -X11.8c >> $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    awk '{if ($4<=10) print $1, $2, $3}' ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean.xyaszs |\
      gmt psxy $proj $rgn -Cage.cpt -St0.4c -O -K >> $ps
    gmt psscale -Dx11.0c/0.5c+w9c/0.25c+ef -Clitho.cpt -B50f10+l"z@-LAB@- (km)" -O -K >> $ps
    echo "0.55 1.05 SL2013sv" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    echo "0.72 -0.05 b" | gmt pstext $rgnx $scalex -F+f18p+jBR -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    gmt psscale -Dx-4.5c/-1.1c+w7c/0.25c+ef+h -Cage.cpt -B5f1+l"Age (Ma)" -O -K >> $ps

    # Plot panel c, mean LAB SL2013
    model="SL2013"
    LABgrid="${fold_input_LAB}/ANT-20/std.grd"
    gmt makecpt -T0/50/2 -Cbroc -D > litho.cpt
    #gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y12.3c -X11.8c >> $ps
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y-12.3c -X-11.8c >> $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    gmt psscale -Dx-1.3c/0.5c+w9c/0.25c+ef+mal -Clitho.cpt -B10f5+l"@~\163@~@-LAB@- (km)" -O -K >> $ps
    echo "0.08 1.05 ANT-20" | gmt pstext $rgnx $scalex -F+f18p+jTR -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    echo "-0.06 -0.05 c" | gmt pstext $rgnx $scalex -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    # Plot panel d, mean LAB CAM2016
    model="CAM2016"
    LABgrid="${fold_input_LAB}/CAM2016/CAM2016Litho.nc"
    gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
    gmt psbasemap $rgn_map $proj_map -Bxa30f15 -K -O -Y0.0c -X11.8c >> $ps
    gmt grdimage $LABgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
    awk '{if ($4<=10) print $1, $2, $3}' ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean.xyaszs |\
      gmt psxy $proj $rgn -Cage.cpt -St0.4c -O -K >> $ps
    gmt psscale -Dx11.0c/0.5c+w9c/0.25c+ef -Clitho.cpt -B50f10+l"z@-LAB@- (km)" -O -K >> $ps
    echo "0.55 1.05 CAM2016" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    echo "0.72 -0.05 d" | gmt pstext $rgnx $scalex -F+f18p+jBR -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    rm litho.cpt

    # Plot panel e, mean LAB ANT-20 vs magmatism age
    model="ANT-20"
    rgn="-R-0.5/23.5/30/100"
    scale="-JX11c/-10c"
    gmt psbasemap $rgn $scale -Bpx5f1+l"Age (Ma)" -Bpy20f10+l"ANT-20 z@-LAB@- (km)" -BWSne -O -K -Y-11.3c -X-12.8c >> $ps
    awk '{print $1, $2, $3, $4}' ${fold_input_conductive_isotherms}/1200C_contour_envelope.txt |\
      gmt psxy $rgn $scale -Gpurple -L+b+t10 -t50 -O -K >> $ps
    awk '{print $1, $2}' ${fold_input_conductive_isotherms}/1200C_contour_envelope.txt |\
      gmt psxy $rgn $scale -W1.0p,71/7/87 -O -K >> $ps
    #awk '{if($4<=10) print $3, $5, $3, $4, $6}' ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean_range.xyaszs | gmt psxy $rgn $scale -E7p/1.0p -St0.4c -O -K -Cage.cpt >> $ps
    awk '{if($4<=10) print $3, $5, $3, $4, $6}' ${fold_data_output}/Antarctica_magmatism_LAB_${model}_blockmean.xyaszs | gmt psxy $rgn $scale -E7p/1.0p,black -St0.4c -O -K -Cage.cpt >> $ps
    echo "0.012 0.057 e" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps

    # Plot panel f, LAB vs magmatism age correlation coefficient
    rgn="-R-0.7/1.0/0/7"
    scale="-JX11c/10c"
    gmt psbasemap $rgn $scale -Bpx0.2f0.1+l"@~\162@~" -Bpy0.5f0.25+l"f(@~\162@~) (%)" -BwSnE -O -K -Y0c -X12.8c >> $ps
    awk '{print $3}' ${fold_data_output}/Antarctica_magmatism_LAB_CAM2016_spearman.spsp | gmt pshistogram -Z1 $rgn $scale -W0.02 -Gred -O -K -t50 >> $ps
    awk '{print $3}' ${fold_data_output}/Antarctica_magmatism_LAB_ANT-20_spearman.spsp | gmt pshistogram -Z1 $rgn $scale -W0.02 -Ggreen4 -O -K -t50 >> $ps
    awk '{print $3}' ${fold_data_output}/Antarctica_magmatism_LAB_SL2013_spearman.spsp | gmt pshistogram -Z1 $rgn $scale -W0.02 -Gblue -O -K -t50 >> $ps
    echo "0.7075 0.068 f" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
    printf "0.296, 0.0\n0.296,8.0" | gmt psxy $rgn $scale -W2.0p,- -O >> $ps

    gmt psconvert -Tj -E600 -A0.1c -P -Z $ps
}

prepare_conductive_cooling_envelope(){

  upper="${fold_input_conductive_isotherms}/subset/1200C_contour_TBL_continent-kradqz-H0.7e-6-s65km-tc10-dt15-1353C-250km.txt"
  lower="${fold_input_conductive_isotherms}/subset/1200C_contour_TBL_continent-kradqz-H0.7e-6-s35km-tc40-dt15-1507C-250km.txt"
  mean="${fold_input_conductive_isotherms}/subset/1200C_contour_TBL_continent-kradqz-H0.7e-6-s50km-tc25-dt15-1468C-250km.txt"
  awk '{print $2}' $lower > junk_lower
  awk '{print $2}' $upper > junk_upper
  paste $mean junk_lower junk_upper > ${fold_input_conductive_isotherms}/1200C_contour_envelope.txt
  rm junk*

}

fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_scripts=${fold_base}/scripts
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')

fold_input_magmatism=$(awk '$1 ~ /^input_magmatism/' config.ini | awk '{print $3}')
fold_input_magmatism=${fold_data_input}/${fold_input_magmatism}
fold_input_LAB=$(awk '$1 ~ /^input_LAB_models/' config.ini | awk '{print $3}')
fold_input_LAB=${fold_data_input}/${fold_input_LAB}
fold_input_conductive_isotherms=$(awk '$1 ~ /^input_conductive_isotherms/' config.ini | awk '{print $3}')
fold_input_conductive_isotherms=${fold_data_input}/${fold_input_conductive_isotherms}

fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')

mkdir -p ${fold_plot_output}
mkdir -p ${fold_data_output}

#prepare_conductive_cooling_envelope
#plot_LAB_magmatism_paper_figure
plot_LAB

rm -f gmt.conf gmt.history *.cpt junk* ints
