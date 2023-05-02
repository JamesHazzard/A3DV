#!/bin/bash

construct_polygonfiles(){

  cat ${fold_data_input_polygons}/ant_full_drainagesystem_polygons.txt > ${fold_data_input_polygons}/polygons.txt
  cat ${fold_data_input_polygons}/extra_polygon.txt >> ${fold_data_input_polygons}/polygons.txt

  no_polygons=$(awk '{if($1=="East") print NF-1}' ${fold_data_input_polygons}/ant_full_drainagesystem_indices.txt)
  i_min=1
  i_max=$no_polygons
  echo $i_min $i_max
  #issue with polygon 1 and polygon 16 (ID 2, 17) due to crossing pole or longitude discontinuity -> split and add polygons
  echo ">" > ${fold_data_input_polygons}/east_ant_polygonfile.xy
  for ((i=$i_min; i<=$i_max; i++)); do

    polygon_ID=$(awk '{if($1=="East") print $('$i'+1)}' ${fold_data_input_polygons}/ant_full_drainagesystem_indices.txt)
    echo constructing polygon ${i} of ${i_max}, polygon ID is $polygon_ID
    awk '{if($3=='$polygon_ID') print $2, $1}' ${fold_data_input_polygons}/polygons.txt | awk '{if($1<0.0) {printf ("%.7f %.7f\n", $1+360.0, $2)}}' > junko1
    l=$(wc -l junko1 | awk '{print $1}')
    awk '{if('$l'>0) print $0}' junko1 >> ${fold_data_input_polygons}/east_ant_polygonfile.xy
    tail -n 1 junko1 | awk '{if('$l'>0) print ">"}' >> ${fold_data_input_polygons}/east_ant_polygonfile.xy
    awk '{if($3=='$polygon_ID') print $2, $1}' ${fold_data_input_polygons}/polygons.txt | awk '{if($1>=0.0) {printf ("%.7f %.7f\n", $1, $2)}}' > junko2
    l=$(wc -l junko2 | awk '{print $1}')
    awk '{if('$l'>0) print $0}' junko2 >> ${fold_data_input_polygons}/east_ant_polygonfile.xy
    tail -n 1 junko2 | awk '{if('$l'>0) print ">"}' >> ${fold_data_input_polygons}/east_ant_polygonfile.xy

  done

  no_polygons=$(awk '{if($1=="West") print NF-1}' ${fold_data_input_polygons}/ant_full_drainagesystem_indices.txt)
  i_min=1
  i_max=$no_polygons
  echo $i_min $i_max
  #issue with polygon 1 and polygon 16 (ID 2, 17) due to crossing pole or longitude discontinuity -> split and add polygons
  echo ">" > ${fold_data_input_polygons}/west_ant_polygonfile.xy
  for ((i=$i_min; i<=$i_max; i++)); do

    polygon_ID=$(awk '{if($1=="West") print $('$i'+1)}' ${fold_data_input_polygons}/ant_full_drainagesystem_indices.txt)
    echo constructing polygon ${i} of ${i_max}, polygon ID is $polygon_ID
    awk '{if($3=='$polygon_ID') print $2, $1}' ${fold_data_input_polygons}/polygons.txt | awk '{if($1<0.0) {printf ("%.7f %.7f\n", $1+360.0, $2)}}' > junko1
    l=$(wc -l junko1 | awk '{print $1}')
    awk '{if('$l'>0) print $0}' junko1 >> ${fold_data_input_polygons}/west_ant_polygonfile.xy
    tail -n 1 junko1 | awk '{if('$l'>0) print ">"}' >> ${fold_data_input_polygons}/west_ant_polygonfile.xy
    awk '{if($3=='$polygon_ID') print $2, $1}' ${fold_data_input_polygons}/polygons.txt | awk '{if($1>=0.0) {printf ("%.7f %.7f\n", $1, $2)}}' > junko2
    l=$(wc -l junko2 | awk '{print $1}')
    awk '{if('$l'>0) print $0}' junko2 >> ${fold_data_input_polygons}/west_ant_polygonfile.xy
    tail -n 1 junko2 | awk '{if('$l'>0) print ">"}' >> ${fold_data_input_polygons}/west_ant_polygonfile.xy
    
  done

  cat ${fold_data_input_polygons}/west_ant_polygonfile.xy > ${fold_data_input_polygons}/ant_polygonfile.xy
  cat ${fold_data_input_polygons}/east_ant_polygonfile.xy >> ${fold_data_input_polygons}/ant_polygonfile.xy

}

divide_grid_ASE(){

  gmt grd2xyz ${file_preamble}.grd | gmt select -F${fold_data_input_polygons}/ASE_polygon.txt > ${file_preamble_out}_ASE.xyz
  gmt xyz2grd ${file_preamble_out}_ASE.xyz -R0/360/-90/-60 -I${inc} -G${file_preamble_out}_ASE.grd
  python3 summary_calculator.py -a ${file_preamble_out}_ASE.xyz

}

divide_grid_Antarctica(){

  gmt grd2xyz ${file_preamble}.grd | gmt select -F${fold_data_input_polygons}/ant_polygonfile.xy > ${file_preamble_out}_ant.xyz
  gmt xyz2grd ${file_preamble_out}_ant.xyz -R0/360/-90/-60 -I${inc} -G${file_preamble_out}_ant.grd
  python3 summary_calculator.py -a ${file_preamble_out}_ant.xyz

  #proj_map="-JS0/-90/10c"
  #rgn_map="-R0/360/-90/-60"
  #gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -M > test_ant_polygonfile.xy
  #gmt grd2xyz ${file_preamble}.grd | gmt select -Ftest_ant_polygonfile.xy > ${file_preamble_out}_ant.xyz
  #gmt xyz2grd ${file_preamble_out}_ant.xyz -R0/360/-90/-60 -I${inc} -G${file_preamble_out}_ant.grd

}

divide_grid_east(){

  if [[ -e ${file_preamble}_east.grd ]]; then
    echo grid already exists
    gmt grd2xyz ${file_preamble}.grd | gmt select -F${fold_data_input_polygons}/east_ant_polygonfile.xy > ${file_preamble_out}_east.xyz
    gmt xyz2grd ${file_preamble_out}_east.xyz -R0/360/-90/-60 -I${inc} -G${file_preamble_out}_east.grd
    rm ${file_preamble_out}_east.xyz
  else
    gmt grd2xyz ${file_preamble}.grd | gmt select -F${fold_data_input_polygons}/east_ant_polygonfile.xy > ${file_preamble_out}_east.xyz
    gmt xyz2grd ${file_preamble_out}_east.xyz -R0/360/-90/-60 -I${inc} -G${file_preamble_out}_east.grd
    rm ${file_preamble_out}_east.xyz
  fi

  gmt grd2xyz ${file_preamble_out}_east.grd > ${file_preamble_out}_east.xyz
  python3 summary_calculator.py -a ${file_preamble_out}_east.xyz

  ps="${fold_plot_output}/maps/${grid_type}_${summary_type}_east_map.ps"
  jpg="${fold_plot_output}/maps/${grid_type}_${summary_type}_east_map.jpg"

  proj="-JS0/-90/10c"
  rgn="-R0/360/-90/-60"

  gmt gmtset COLOR_NAN gray PS_MEDIA a0 MAP_FRAME_TYPE plain
  grid_full=${file_preamble_out}_east.grd
  echo $grid_type ${grid_type}
  if [[ "${grid_type}" == "HF" ]]; then
    echo $grid_type HF
    gmt makecpt -T40/100/2.5 -Cvik -D > litho.cpt
  elif [[ "${grid_type}" == "LAB1200" ]]; then
    echo $grid_type LAB1200
    gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
  elif [[ "${grid_type}" == "eta" ]]; then
    echo $grid_type eta
    gmt makecpt -Chot -D -T18/23/0.25 -I -G0.08/1 > litho.cpt
  fi
  gmt grdimage $grid_full $proj $rgn -Clitho.cpt -B30f10g30 -K -E600 -X12.5c -Y25c > $ps
  gmt pscoast $proj $rgn -Dh -W0.5p -A50000/0/0 -O -K >> $ps
  if [[ "${grid_type}" == "HF" ]]; then
    gmt psscale -D+5.125c/-1.5c+w10c/0.2c+jMC+h+e -Clitho.cpt -B10f10+l"@~\155@~@-HF@- (mW m@+-2@+)" -O >> $ps
  elif [[ "${grid_type}" == "LAB1200" ]]; then
    gmt psscale -D+5.125c/-1.5c+w10c/0.2c+jMC+h+e -Clitho.cpt -B50f10+l"z@-LAB@- (km)" -O >> $ps
  elif [[ "${grid_type}" == "eta" ]]; then
    gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Clitho.cpt -B1f0.25+l"log@-10@-@~\150@~ (Pa s)" -O >> $ps
  fi

  gmt psconvert $ps -Tj -E600 -A0.1c/0.1c -P -Z

  rm gmt.* litho.cpt
}

divide_grid_west(){

  if [[ -e ${file_preamble}_west.grd ]]; then
    echo grid already exists
    gmt grd2xyz ${file_preamble}.grd | gmt select -F${fold_data_input_polygons}/west_ant_polygonfile.xy > ${file_preamble_out}_west.xyz
    gmt xyz2grd ${file_preamble_out}_west.xyz -R0/360/-90/-60 -I${inc} -G${file_preamble_out}_west.grd
    rm ${file_preamble_out}_west.xyz
  else
    gmt grd2xyz ${file_preamble}.grd | gmt select -F${fold_data_input_polygons}/west_ant_polygonfile.xy > ${file_preamble_out}_west.xyz
    gmt xyz2grd ${file_preamble_out}_west.xyz -R0/360/-90/-60 -I${inc} -G${file_preamble_out}_west.grd
    rm ${file_preamble_out}_west.xyz
  fi

  gmt grd2xyz ${file_preamble_out}_west.grd > ${file_preamble_out}_west.xyz
  python3 summary_calculator.py -a ${file_preamble_out}_west.xyz

  ps="${fold_plot_output}/maps/${grid_type}_${summary_type}_west_map.ps"
  jpg="${fold_plot_output}/maps/${grid_type}_${summary_type}_west_map.jpg"

  proj="-JS0/-90/10c"
  rgn="-R0/360/-90/-60"

  gmt gmtset COLOR_NAN gray PS_MEDIA a0 MAP_FRAME_TYPE plain
  grid_full=${file_preamble_out}_west.grd
  echo $grid_type ${grid_type}
  if [[ "${grid_type}" == "HF" ]]; then
    echo $grid_type HF
    gmt makecpt -T40/100/2.5 -Cvik -D > litho.cpt
  elif [[ "${grid_type}" == "LAB1200" ]]; then
    echo $grid_type LAB1200
    gmt makecpt -T0/350/10 -Cbroc -D > litho.cpt
  elif [[ "${grid_type}" == "eta" ]]; then
    echo $grid_type eta
    gmt makecpt -Chot -D -T18/23/0.25 -I -G0.08/1 > litho.cpt
  fi
  gmt grdimage $grid_full $proj $rgn -Clitho.cpt -B30f10g30 -K -E600 -X12.5c -Y25c > $ps
  gmt pscoast $proj $rgn -Dh -W0.5p -A50000/0/0 -O -K >> $ps
  if [[ "${grid_type}" == "HF" ]]; then
    gmt psscale -D+5.125c/-1.5c+w10c/0.2c+jMC+h+e -Clitho.cpt -B10f10+l"@~\155@~@-HF@- (mW m@+-2@+)" -O >> $ps
  elif [[ "${grid_type}" == "LAB1200" ]]; then
    gmt psscale -D+5.125c/-1.5c+w10c/0.2c+jMC+h+e -Clitho.cpt -B50f10+l"z@-LAB@- (km)" -O >> $ps
  elif [[ "${grid_type}" == "eta" ]]; then
    gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Clitho.cpt -B1f0.25+l"log@-10@-@~\150@~ (Pa s)" -O >> $ps
  fi

  gmt psconvert $ps -Tj -E600 -A0.1c/0.1c -P -Z

  rm gmt.* litho.cpt
}

plot_LAB_histogram_panel(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  ps="${fold_plot_output}/histograms/LAB_histogram_panel.ps"
  jpg="${fold_plot_output}/histograms/LAB_histogram_panel.jpg"

  max_count=2500
  rgn="-R30/350/0/${max_count}"
  scale="-JX10c/10c"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  grid_type="LAB1200"

  gmt psbasemap $rgn $scale -Bpx50f10+l"z@-LAB@- (km)" -Bpy500f100+l"f(z@-LAB@-)" -BwSnE -K -Y25c -X25c > $ps
  gmt grd2xyz ${fold_data_output}/${grid_type}_mean_filtered_east.grd | awk '{print $3}' | grep -v "NaN" | gmt pshistogram -Z0 $rgn $scale -W10 -G28/99/114 -O -K -t30 >> $ps
  gmt grd2xyz ${fold_data_output}/${grid_type}_mean_filtered_west.grd | awk '{print $3}' | grep -v "NaN" | gmt pshistogram -Z0 $rgn $scale -W10 -G204/160/24 -O -K -t30 >> $ps
  median=$(awk '{if(NR==3) print $1}' ${fold_data_output}/${grid_type}_${summary_type}_filtered_west.xyz.summary)
  MAD=$(awk '{if(NR==4) print $1}' ${fold_data_output}/${grid_type}_${summary_type}_filtered_west.xyz.summary)
  echo $median $MAD | awk '{print $1-$2, $1, $1+$2}'
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1-$2, 0.0, $1-$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,163/99/33,dashed -O -K >> $ps
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1+$2, 0.0, $1+$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,163/99/33,dashed -O -K >> $ps
  median=$(awk '{if(NR==3) print $1}' ${fold_data_output}/${grid_type}_${summary_type}_filtered_east.xyz.summary)
  MAD=$(awk '{if(NR==4) print $1}' ${fold_data_output}/${grid_type}_${summary_type}_filtered_east.xyz.summary)
  echo $median $MAD | awk '{print $1-$2, $1, $1+$2}'
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1-$2, 0.0, $1-$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,13/63/74,dashed -O -K >> $ps
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1+$2, 0.0, $1+$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,13/63/74,dashed -O -K >> $ps
  awk '{if(NR==3) printf "%.7f %.7f\n%.7f %.7f", $1, 0.0, $1, '${max_count}'}' ${file_preamble_out}_west.xyz.summary | gmt psxy $rgn $scale -W1.5p,163/99/33,dashed -O -K >> $ps
  awk '{if(NR==3) printf "%.7f %.7f\n%.7f %.7f", $1, 0.0, $1, '${max_count}'}' ${file_preamble_out}_east.xyz.summary | gmt psxy $rgn $scale -W1.5p,13/63/74,dashed -O -K >> $ps
  echo "0.6525 0.9775 b" | gmt pstext $rgnx $scalex -F+f18p+jTR -C+tO -Gwhite -W1p,black -N -O -K >> $ps

  proj="-JS0/-90/10c"
  rgn="-R0/360/-90/-61"

  gmt psbasemap $rgn $proj -Bxa30f15 -O -K -Y0c -X-11c >> $ps
  gmt pscoast $rgn $proj -Dh -A50000/0/2 -G28/99/114 -W0p -O -K >> $ps
  #gmt pscoast $rgn $proj -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psxy ${fold_data_input_polygons}/west_ant_polygonfile.xy $rgn $proj -G204/160/24 -O -K -L >> $ps
  count_systems=$(tail -n 1 ${fold_data_input_polygons}/ant_full_drainagesystem_polygons.txt | awk '{print $3}')
  #count_systems=1
  for ((i=1; i<=$count_systems; i++)); do
    awk '{if($3=='$i') print $2, $1}' ${fold_data_input_polygons}/ant_full_drainagesystem_polygons.txt | gmt psxy $rgn $proj -W0.5p,black -L -O -K >> $ps
    echo $i
  done
  gmt pscoast $rgn $proj -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  echo "-0.05 0.9775 a" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O >> $ps

  gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

plot_HF(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  # Set plot name
  ps="${fold_plot_output}/maps/HF.ps"
  jpg="${fold_plot_output}/maps/HF.jpg"

  # Set plot area
  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  grid_type="HF"

  # Plot panel a, mean HF
  HFgrid=${fold_data_output}/${grid_type}_mean_filtered_ant.grd
  gmt makecpt -T30/110/2.5 -Cvik -D > litho.cpt
  gmt psbasemap $rgn_map $proj_map -B0 -K -Y25c -X25c > $ps
  gmt grdimage $HFgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Clitho.cpt -B10f10+l"@~\155@~@-HF@- (mW m@+-2@+)" -Al -O -K >> $ps
  echo "0.03 1.05 a" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K >> $ps
  rm litho.cpt

  # Plot panel b, std HF
  HFgrid=${fold_data_output}/${grid_type}_std_filtered_ant.grd
  gmt makecpt -T0/10/0.5 -Cvik -D > litho.cpt
  gmt psbasemap $rgn_map $proj_map -B0 -K -O -Y0c -X12c >> $ps
  gmt grdimage $HFgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Clitho.cpt -B2f1+l"@~\163@~@-HF@- (mW m@+-2@+)" -Al -O -K >> $ps
  echo "0.63 1.05 b" | gmt pstext $rgnx $scalex -F+f18p+jTR -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O >> $ps
  rm litho.cpt

  gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

plot_HF_histogram_panel(){

  # Set GMT plotting parameters
  gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

  # Set plot name
  ps="${fold_plot_output}/histograms/HF_histogram_panel.ps"
  jpg="${fold_plot_output}/histograms/HF_histogram_panel.jpg"

  # Set plot area
  proj_map="-JS0/-90/10c"
  rgn_map="-R0/360/-90/-60"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  grid_type="HF"

  # Plot panel a, mean HF
  HFgrid=${fold_data_output}/${grid_type}_mean_filtered_ant.grd
  gmt makecpt -T30/110/2.5 -Cvik -D > litho.cpt
  gmt psbasemap $rgn_map $proj_map -B0 -K -Y25c -X25c > $ps
  gmt grdimage $HFgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Clitho.cpt -B10f10+l"@~\155@~@-HF@- (mW m@+-2@+)" -Al -O -K >> $ps
  echo "0.03 1.05 a" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K >> $ps
  rm litho.cpt

  # Plot panel b, std HF
  HFgrid=${fold_data_output}/${grid_type}_std_filtered_ant.grd
  gmt makecpt -T0/10/0.5 -Cvik -D > litho.cpt
  gmt psbasemap $rgn_map $proj_map -B0 -K -O -Y0c -X12c >> $ps
  gmt grdimage $HFgrid -Clitho.cpt $rgn_map $proj_map -E600 -O -K >> $ps
  gmt pscoast $rgn_map $proj_map -Dh -A50000/0/2 -Wthin,black -O -K >> $ps
  gmt psscale -Dx0.5c/-2.0c+w9c/0.25c+e+h -Clitho.cpt -B2f1+l"@~\163@~@-HF@- (mW m@+-2@+)" -Al -O -K >> $ps
  echo "0.63 1.05 b" | gmt pstext $rgnx $scalex -F+f18p+jTR -C+tO -Gwhite -W1p,black -N -O -K >> $ps
  gmt psbasemap $rgn_map $proj_map -Bxa30f15 -O -K >> $ps
  rm litho.cpt

  max_count=920
  rgn="-R30/110/0/${max_count}"
  scale="-JX10c/10c"
  rgnx="-R0/1/0/1"
  scalex="-JX15c/10c"

  gmt psbasemap $rgn $scale -Bpx10f10+l"@~\155@~@-HF@- (mW m@+-2@+)" -Bpy100f50+l"f(@~\155@~@-HF@-)" -BwSnE -O -K -Y-13.2c -X-0.5c >> $ps
  gmt grd2xyz ${fold_data_output}/${grid_type}_mean_filtered_east.grd | awk '{print $3}' | grep -v "NaN" | gmt pshistogram -Z0 $rgn $scale -W2.5 -G28/99/114 -O -K -t30 >> $ps
  gmt grd2xyz ${fold_data_output}/${grid_type}_mean_filtered_west.grd | awk '{print $3}' | grep -v "NaN" | gmt pshistogram -Z0 $rgn $scale -W2.5 -G204/160/24 -O -K -t30 >> $ps
  median=$(awk '{if(NR==3) print $1}' ${fold_data_output}/${grid_type}_mean_filtered_west.xyz.summary)
  MAD=$(awk '{if(NR==4) print $1}' ${fold_data_output}/${grid_type}_mean_filtered_west.xyz.summary)
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1-$2, 0.0, $1-$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,163/99/33,dashed -O -K >> $ps
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1+$2, 0.0, $1+$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,163/99/33,dashed -O -K >> $ps
  median=$(awk '{if(NR==3) print $1}' ${fold_data_output}/${grid_type}_mean_filtered_east.xyz.summary)
  MAD=$(awk '{if(NR==4) print $1}' ${fold_data_output}/${grid_type}_mean_filtered_east.xyz.summary)
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1-$2, 0.0, $1-$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,13/63/74,dashed -O -K >> $ps
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1+$2, 0.0, $1+$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,13/63/74,dashed -O -K >> $ps
  awk '{if(NR==3) printf "%.7f %.7f\n%.7f %.7f", $1, 0.0, $1, '${max_count}'}' ${fold_data_output}/${grid_type}_mean_filtered_west.xyz.summary | gmt psxy $rgn $scale -W1.5p,163/99/33,dashed -O -K >> $ps
  awk '{if(NR==3) printf "%.7f %.7f\n%.7f %.7f", $1, 0.0, $1, '${max_count}'}' ${fold_data_output}/${grid_type}_mean_filtered_east.xyz.summary | gmt psxy $rgn $scale -W1.5p,13/63/74,dashed -O -K >> $ps
  echo "0.6525 0.9775 d" | gmt pstext $rgnx $scalex -F+f18p+jTR -C+tO -Gwhite -W1p,black -N -O -K >> $ps

  max_count=8000
  rgn="-R30/110/0/${max_count}"
  scale="-JX10c/10c"

  gmt psbasemap $rgn $scale -Bpx10f10+l"@~\155@~@-HF@- (mW m@+-2@+)" -Bpy1000f500+l"f(@~\155@~@-HF@-)" -BWSne -O -K -Y0c -X-11.0c >> $ps
  gmt grd2xyz ${fold_data_output}/${grid_type}_mean_filtered_east.grd | awk '{print $3}' | grep -v "NaN" | gmt pshistogram -Z0 $rgn $scale -W2.5 -G28/99/114 -O -K -t30 >> $ps
  gmt grd2xyz ${fold_data_output}/${grid_type}_mean_filtered_west.grd | awk '{print $3}' | grep -v "NaN" | gmt pshistogram -Z0 $rgn $scale -W2.5 -G204/160/24 -O -K -t30 >> $ps
  median=$(awk '{if(NR==3) print $1}' ${fold_data_output}/${grid_type}_mean_filtered_west.xyz.summary)
  MAD=$(awk '{if(NR==4) print $1}' ${fold_data_output}/${grid_type}_mean_filtered_west.xyz.summary)
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1-$2, 0.0, $1-$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,163/99/33,dashed -O -K >> $ps
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1+$2, 0.0, $1+$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,163/99/33,dashed -O -K >> $ps
  median=$(awk '{if(NR==3) print $1}' ${fold_data_output}/${grid_type}_mean_filtered_east.xyz.summary)
  MAD=$(awk '{if(NR==4) print $1}' ${fold_data_output}/${grid_type}_mean_filtered_east.xyz.summary)
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1-$2, 0.0, $1-$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,13/63/74,dashed -O -K >> $ps
  echo $median $MAD | awk '{printf "%.7f %.7f\n%.7f %.7f", $1+$2, 0.0, $1+$2, '${max_count}'}' | gmt psxy $rgn $scale -W1.0p,13/63/74,dashed -O -K >> $ps
  awk '{if(NR==3) printf "%.7f %.7f\n%.7f %.7f", $1, 0.0, $1, '${max_count}'}' ${fold_data_output}/${grid_type}_mean_filtered_west.xyz.summary | gmt psxy $rgn $scale -W1.5p,163/99/33,dashed -O -K >> $ps
  awk '{if(NR==3) printf "%.7f %.7f\n%.7f %.7f", $1, 0.0, $1, '${max_count}'}' ${fold_data_output}/${grid_type}_mean_filtered_east.xyz.summary | gmt psxy $rgn $scale -W1.5p,13/63/74,dashed -O -K >> $ps
  echo "0.0125 0.9775 c" | gmt pstext $rgnx $scalex -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O >> $ps
  gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

divide_HF(){
depth=""
grid_type="HF"
summary_type="mean"
if [[ ${depth} == "" ]]; then
  grid=${grid_type}
else
  grid=${depth}_km_${grid_type}
fi
echo $grid
divide_grid_west
divide_grid_east
divide_grid_ASE
divide_grid_Antarctica
summary_type="std"
if [[ ${depth} == "" ]]; then
  grid=${grid_type}
else
  grid=${depth}_km_${grid_type}
fi
echo $grid
divide_grid_west
divide_grid_east
divide_grid_ASE
divide_grid_Antarctica
}

divide_LAB(){
depth=""
grid_type="LAB"
summary_type="mean"
if [[ ${depth} == "" ]]; then
  grid=${grid_type}
else
  grid=${depth}_km_${grid_type}
fi
echo $grid
#divide_grid_west
#divide_grid_east
divide_grid_ASE
summary_type="std"
if [[ ${depth} == "" ]]; then
  grid=${grid_type}
else
  grid=${depth}_km_${grid_type}
fi
echo $grid
#divide_grid_west
#divide_grid_east
divide_grid_ASE
}

fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')
Tp=$(awk '$1 ~ /^potential_temperature/' config.ini | awk '{print $3}')
sol50=$(awk '$1 ~ /^solidus_50km/' config.ini | awk '{print $3}')
fold_data_output=${fold_data_output}/Tp_${Tp}_sol50_${sol50}
fold_plot_output=${fold_plot_output}/Tp_${Tp}_sol50_${sol50}
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
input_polygons=$(awk '$1 ~ /^input_polygons/' config.ini | awk '{print $3}')
input_geotherm_grids=$(awk '$1 ~ /^input_geotherm_grids/' config.ini | awk '{print $3}')
fold_data_input_polygons=${fold_data_input}/${input_polygons}
fold_data_input_geotherm_grids=${fold_data_input}/${input_geotherm_grids}

mkdir -p ${fold_data_output}
mkdir -p ${fold_plot_output}
mkdir -p ${fold_plot_output}/maps
mkdir -p ${fold_plot_output}/histograms

inc="0.5d"

#construct_polygonfiles

grid_type="LAB1200"
summary_type="mean"
file_preamble=${fold_data_input_geotherm_grids}/Tp_${Tp}_sol50_${sol50}/${grid_type}/distribution/summary/${grid_type}_${summary_type}_filtered
file_preamble_out=${fold_data_output}/${grid_type}_${summary_type}_filtered
#divide_grid_west
#divide_grid_east
#divide_grid_ASE
#divide_grid_Antarctica

#plot_LAB_histogram_panel

grid_type="HF"
summary_type="mean"
file_preamble=${fold_data_input_geotherm_grids}/Tp_${Tp}_sol50_${sol50}/${grid_type}/distribution/summary/${grid_type}_${summary_type}_filtered
file_preamble_out=${fold_data_output}/${grid_type}_${summary_type}_filtered
#divide_grid_west
#divide_grid_east
#divide_grid_ASE
#divide_grid_Antarctica
summary_type="std"
file_preamble=${fold_data_input_geotherm_grids}/Tp_${Tp}_sol50_${sol50}/${grid_type}/distribution/summary/${grid_type}_${summary_type}_filtered
file_preamble_out=${fold_data_output}/${grid_type}_${summary_type}_filtered
#divide_grid_west
#divide_grid_east
#divide_grid_ASE
#divide_grid_Antarctica

#plot_HF_histogram_panel
plot_HF
