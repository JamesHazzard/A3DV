#!/bin/bash

plot_1d_viscosity_profile(){

    loc="ASE"

    # Set GMT plotting parameters   
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

    plot_name="${fold_plot_output}/viscosity_1d_profile_${loc}"
    ps="${plot_name}.ps"
    jpg="${plot_name}.jpg"

    rgn="-R18/24/75/325"
    scale="-JX7.5c/-15c"

    gmt psbasemap $rgn $scale -B0 -X20c -Y20c -K > $ps
    # create the viscosity geotherms
    # sample region over finite radius
    # average viscosity over lateral spatial domain
    # plot viscosity as a function of depth

    for ((idx=0; idx<10; idx++)); do

        echo "calculating 1d viscosity profile at location ${loc} for anelasticity model ${idx}.."

        f_out=${input_viscosity}/location_${loc}/1d_profile_Tp_${Tp}_sol50_${sol50}_idx_${idx}_${loc}.txt
        if [[ -e ${f_out} ]]; then
            echo "exists, skipping"
        else
            echo -n > ${f_out}
            for z in $(seq 100 25 300); do

                f_viscosity=${input_viscosity}/location_${loc}/depth_${z}_Tp_${Tp}_sol50_${sol50}_idx_${idx}_${loc}.grd
                gmt grdmath $f_viscosity MEAN = junk
                gmt grd2xyz junk > junk2
                visc=$(head -n 1 junk2 | awk '{print $3}')
                echo $z $visc >> ${f_out}

            done
        fi

        awk '{if($1<=300) print $2, $1}' $f_out | gmt psxy $rgn $scale -W1p,black,solid -O -K >> $ps

    done

    gmt psbasemap $rgn $scale -Bx1f0.5+l"log@-10@-@~\150@~" -By50f25+l"Depth (km)" -BWSne -O >> $ps
    gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

    rm -f gmt.*

}

plot_kde_presentation(){

    # set reasonable estimates of geodetic 'observed' viscosity (fitted to observed uplift rates - sig. trade-off with LAB depth)
    echo -n > junk
    echo "ASE B18 18.4 19.4 150 175 150 250 c 30c 30c NesW Nesw 0 3.5" >> junk
    echo "AP S20 17.5 19.0 125 175 125 250 d 6c 0c Nesw NEsw 0 6.5" >> junk
    echo "AP I11 19.3 20.0 125 175 125 250 e -16c -9c neSW neSw 0 2.5" >> junk
    echo "AP Wo15 20.0 20.5 125 175 125 250 f 6c 0c neSw nESw 0 2.5" >> junk

    # Set GMT plotting parameters   
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain FORMAT_GEO_OUT F FORMAT_GEO_MAP F MAP_ANNOT_OBLIQUE 15

    fold_plot_output_histogram=${fold_plot_output}
    mkdir -p ${fold_plot_output_histogram}
    plot_name="${fold_plot_output_histogram}/viscosity_kde_ss_and_band"
    ps="${plot_name}.ps"
    jpg="${plot_name}.jpg"

    rgnx="-R0/1/0/1"
    scale="-JX10c/8c"
    scale_sub="-JX5c/8c"

    rgn_map="-R-84/-65/-36/-70r"
    proj_map="-JS0/-90/10c"

    idx_line=$(echo 1)
    idx_line_max=$(cat junk | wc -l)

    for ((idx_i=0; idx_i<=1; idx_i++)); do

        for ((idx_j=0; idx_j<=1; idx_j++)); do

            #move_x=$(echo $idx_i $idx_j | awk '{if($2>0) print "11c"; else if ($1<1 && $2<1) print "0c"; else print "-11c"}')
            #move_y=$(echo $idx_i $idx_j | awk '{if($1>0 && $2<1) print "-11c"; else print "0c"}')
            move_x=$(awk '{if(NR=='${idx_line}') print $10}' junk)
            move_y=$(awk '{if(NR=='${idx_line}') print $11}' junk)
            axis_a=$(awk '{if(NR=='${idx_line}') print $12}' junk)
            axis_b=$(awk '{if(NR=='${idx_line}') print $13}' junk)
            axis_x=$(echo $idx_i $idx_j | awk '{print "ns"}')
            axis_y=$(echo $idx_i $idx_j | awk '{if($2<1) print "We"; else print "wE"}')

            # set key variables incl. location, loading model, geodetic viscosity estimate, and file name for apparent viscosity data
            loc=$(awk '{if(NR=='${idx_line}') print $1}' junk)
            loading_model=$(awk '{if(NR=='${idx_line}') print $2}' junk)
            echo $loc $loading_model $idx_i $idx_j
            eta_min=$(awk '{if(NR=='${idx_line}') print $3}' junk)
            eta_max=$(awk '{if(NR=='${idx_line}') print $4}' junk)

            depthtop_a=$(awk '{if(NR=='${idx_line}') print $5}' junk)
            depthbottom_a=$(awk '{if(NR=='${idx_line}') print $6}' junk)
            data_a=${fold_data_output}/apparent_viscosity/distribution/summary/apparent_viscosity_location_${loc}_loading_${loading_model}_depthtop_${depthtop_a}_depthbottom_${depthbottom_a}_kde.txt
            depthtop_b=$(awk '{if(NR=='${idx_line}') print $7}' junk)
            depthbottom_b=$(awk '{if(NR=='${idx_line}') print $8}' junk)
            data_b=${fold_data_output}/apparent_viscosity/distribution/summary/apparent_viscosity_location_${loc}_loading_${loading_model}_depthtop_${depthtop_b}_depthbottom_${depthbottom_b}_kde.txt
            label=$(awk '{if(NR=='${idx_line}') print $9}' junk)

            # set histogram plot region and label
            y_min=17.5
            y_max=21.0
            y_major_tick=1
            y_minor_tick=0.1
            bin_width=0.05
            #x_min=0
            #x_max=$(awk 'BEGIN{a=0}{if($3>0+a) a=$3} END{print 1.2*a}' $data_a)
            x_min=$(awk '{if(NR=='${idx_line}') print $14}' junk)
            x_max=$(awk '{if(NR=='${idx_line}') print $15}' junk)
            x_major_tick=1
            x_minor_tick=0.5
            x_label="f(log@-10@-@~\150@~)"
            y_label="log@-10@-@~\150@~"
            rgn="-R${x_min}/${x_max}/${y_min}/${y_max}"
            
            # plot histogram
            gmt_flag=$(echo $idx_line | awk '{if($1==1) print "-K"; else print "-O -K"}')
            echo $idx_line $gmt_flag
            gmt psbasemap $rgn $scale_sub -Bpy${y_major_tick}f${y_minor_tick}+l"${y_label}" -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -X$move_x -Y$move_y -B${axis_a} $gmt_flag >> $ps
            echo 0 ${eta_min} > junk2
            echo 0 ${eta_max} >> junk2
            echo ${x_max} ${eta_max} >> junk2
            echo ${x_max} ${eta_min} >> junk2
            head -n 1 junk2 >> junk2
            awk '{print $1, $2}' junk2 | gmt psxy $rgn $scale -Ggray -t50 -O -K >> $ps
            dx=9.85
            dy=-7.8
            echo "0 1 (${loc}, ${loading_model})" | gmt pstext $rgnx $scale -F+jBR -D${dx}/${dy} -O -K >> $ps
            dx=0.15
            dy=-0.2
            echo "0 1 ${depthtop_a}-${depthbottom_a} km" | gmt pstext $rgnx $scale -F+jTL -D${dx}/${dy} -O -K >> $ps
            awk '{print $2, $1}' $data_a | gmt psxy $rgn $scale_sub -W1.5p,blue,solid -O -K -t50 >> $ps
            #awk '{print $3, $1}' $data_a | gmt psxy $rgn $scale_sub -W1.5p,red,solid -O -K -t50 >> $ps

            #x_min=0
            #x_max=$(awk 'BEGIN{a=0}{if($3>0+a) a=$3} END{print 1.5*a}' $data_b)
            #y_label="log@-10@-@~\150@~"
            #rgn="-R${x_min}/${x_max}/${y_min}/${y_max}"
            gmt psbasemap $rgn $scale_sub -Bpy${y_major_tick}f${y_minor_tick}+l"${y_label}" -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -X5c -Y0c -B${axis_b} -O -K >> $ps
            awk '{print $2, $1}' $data_b | gmt psxy $rgn $scale_sub -W1.5p,blue,solid -O -K -t50 >> $ps
            #awk '{print $3, $1}' $data_b | gmt psxy $rgn $scale_sub -W1.5p,red,solid -O -K -t50 >> $ps
            dx=4.85
            dy=-0.2
            echo "0 1 ${depthtop_b}-${depthbottom_b} km" | gmt pstext $rgnx $scale -F+jTR -D${dx}/${dy} -O -K >> $ps

            idx_line=$(echo $idx_line | awk '{print $1+1}')

        done

    done

    gmt psbasemap $rgn $scale_sub -B0 -O >> $ps
    gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

plot_geodetic_vs_model_viscosity_kdes_ellipses(){

    # set reasonable estimates of geodetic 'observed' viscosity (fitted to observed uplift rates - sig. trade-off with LAB depth)
    echo -n > junk
    echo "ASE B18 18.4 19.4 150 175 150 250 c -0c -12.5c NesW Nesw 0 3.5" >> junk
    echo "S21AP_perimeter S20 17.5 19.0 125 175 125 250 d 6c 0c Nesw NEsw 0 6.5" >> junk
    echo "I11AP_perimeter I11 19.3 20.0 125 175 125 250 e -16c -9c neSW neSw 0 2.5" >> junk
    echo "Wo15AP_perimeter Wo15 20.0 20.5 125 175 125 250 f 6c 0c neSw nESw 0 2.5" >> junk

    # Set GMT plotting parameters   
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain FORMAT_GEO_OUT F FORMAT_GEO_MAP F MAP_ANNOT_OBLIQUE 15

    fold_plot_output_histogram=${fold_plot_output}
    mkdir -p ${fold_plot_output_histogram}
    plot_name="${fold_plot_output_histogram}/viscosity_kdes_GPS_Ju_spatial_var_ellipses"
    ps="${plot_name}.ps"
    jpg="${plot_name}.jpg"

    rgnx="-R0/1/0/1"
    scale="-JX10c/8c"
    scale_sub="-JX5c/8c"

    rgn_map="-R-84/-65/-36/-70r"
    proj_map="-JS0/-90/10c"
    gmt psbasemap $rgn_map $proj_map -X30c -Y30c -B0 -K > $ps
    gmt makecpt -Chot -D -T18/23/0.1 -I -G0.08/1 > visc.cpt
    logetagrid=${fold_data_input_thermo}/${date}/Tp_${Tp}_sol50_${sol50}/viscosity_grids/distribution/summary/viscosity_depth_150_Tp_${Tp}_sol50_${sol50}_mean.grd
    gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A5000/0/4 -W1.2pt,black -O -K >> $ps
    #awk '{if(NR>=4) print $2, $3, "100k"}' ${fold_data_input_locs}/sample_locs.txt |\
    # gmt psxy $rgn_map $proj_map -SE- -W2p,blue,solid -O -K >> $ps
    for plot_loc in I11 S20 Wo15; do
        awk '{if(NR>1) print $1, $2}' ${fold_data_input_locs}/${plot_loc}_perimeter.txt > junk_perimeter
        head -n 1 junk_perimeter >> junk_perimeter
        gmt psxy junk_perimeter $rgn_map $proj_map -W2p,blue,solid -O -K >> $ps
    done
    dx=9.85
    dy=-9.8
    #gmt psbasemap $rgn_map $proj_map -Bxa5f2.5 -Bpya2f1 -BNESw -O -K >> $ps
    gmt psbasemap $rgn_map $proj_map -Bx15f1g15 -By5f1g5 -BNESw -O -K >> $ps
    echo "0 1 b" | gmt pstext $rgnx $proj_map -F+jBR -W1p,black,solid -Gwhite -TO -D${dx}/${dy} -O -K >> $ps
    gmt psscale -Dx-4.5c/-2.0c+w8c/0.25c+e+h -Cvisc.cpt -B1f0.5+l"log@-10@-@~\155@~@-@~\150@~@- (Pa s)" -Al -O -K >> $ps

    rgn_map="-R-117/-66/-88/-79.5r"
    proj_map="-JS0/-90/10c"
    gmt psbasemap $rgn_map $proj_map -X-11c -B0 -O -K >> $ps
    gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A5000/0/4 -W1.2pt,black -O -K >> $ps
    awk '{if(NR==2) print $2, $3, "100k"}' ${fold_data_input_locs}/sample_locs.txt |\
     gmt psxy $rgn_map $proj_map -SE- -W2p,blue,solid -O -K >> $ps
    dx=0.15
    dy=-9.8
    gmt psbasemap $rgn_map $proj_map -Bx15f1g15 -By5f1g5 -BNeSW -O -K >> $ps
    echo "0 1 a" | gmt pstext $rgnx $proj_map -F+jBL -W1p,black,solid -Gwhite -TO -D${dx}/${dy} -O -K >> $ps

    idx_line=$(echo 1)
    idx_line_max=$(cat junk | wc -l)

    for ((idx_i=0; idx_i<=1; idx_i++)); do

        for ((idx_j=0; idx_j<=1; idx_j++)); do

            move_x=$(awk '{if(NR=='${idx_line}') print $10}' junk)
            move_y=$(awk '{if(NR=='${idx_line}') print $11}' junk)
            axis_a=$(awk '{if(NR=='${idx_line}') print $12}' junk)
            axis_b=$(awk '{if(NR=='${idx_line}') print $13}' junk)
            axis_x=$(echo $idx_i $idx_j | awk '{print "ns"}')
            axis_y=$(echo $idx_i $idx_j | awk '{if($2<1) print "We"; else print "wE"}')

            # set key variables incl. location, loading model, geodetic viscosity estimate, and file name for apparent viscosity data
            loc=$(awk '{if(NR=='${idx_line}') print $1}' junk)
            loading_model=$(awk '{if(NR=='${idx_line}') print $2}' junk)
            echo $loc $loading_model $idx_i $idx_j
            eta_min=$(awk '{if(NR=='${idx_line}') print $3}' junk)
            eta_max=$(awk '{if(NR=='${idx_line}') print $4}' junk)

            depthtop_a=$(awk '{if(NR=='${idx_line}') print $5}' junk)
            depthbottom_a=$(awk '{if(NR=='${idx_line}') print $6}' junk)
            data_a=${fold_data_output}/apparent_viscosity/distribution/summary/apparent_viscosity_location_${loc}_loading_${loading_model}_depthtop_${depthtop_a}_depthbottom_${depthbottom_a}_kde.txt
            depthtop_b=$(awk '{if(NR=='${idx_line}') print $7}' junk)
            depthbottom_b=$(awk '{if(NR=='${idx_line}') print $8}' junk)
            data_b=${fold_data_output}/apparent_viscosity/distribution/summary/apparent_viscosity_location_${loc}_loading_${loading_model}_depthtop_${depthtop_b}_depthbottom_${depthbottom_b}_kde.txt
            label=$(awk '{if(NR=='${idx_line}') print $9}' junk)

            # set histogram plot region and label
            y_min=17.5
            y_max=21.0
            y_major_tick=1
            y_minor_tick=0.1
            bin_width=0.05
            x_min=$(awk '{if(NR=='${idx_line}') print $14}' junk)
            x_max=$(awk '{if(NR=='${idx_line}') print $15}' junk)
            x_major_tick=1
            x_minor_tick=0.5
            x_label="f(log@-10@-@~\150@~)"
            y_label="log@-10@-@~\150@~"
            rgn="-R${x_min}/${x_max}/${y_min}/${y_max}"
            
            # plot histogram
            gmt psbasemap $rgn $scale_sub -Bpy${y_major_tick}f${y_minor_tick}+l"${y_label}" -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -X$move_x -Y$move_y -B${axis_a} -O -K >> $ps
            echo 0 ${eta_min} > junk2
            echo 0 ${eta_max} >> junk2
            echo ${x_max} ${eta_max} >> junk2
            echo ${x_max} ${eta_min} >> junk2
            head -n 1 junk2 >> junk2
            awk '{print $1, $2}' junk2 | gmt psxy $rgn $scale -Ggray -t50 -O -K >> $ps
            dx=9.85
            dy=-7.8
            echo "0 1 (${loc}, ${loading_model})" | gmt pstext $rgnx $scale -F+jBR -D${dx}/${dy} -O -K >> $ps
            dx=0.15
            dy=-7.8
            echo "0 1 ${label}" | gmt pstext $rgnx $scale -F+jBL -W1p,black,solid -Gwhite -TO -D${dx}/${dy} -O -K >> $ps
            dx=0.15
            dy=-0.2
            echo "0 1 ${depthtop_a}-${depthbottom_a} km" | gmt pstext $rgnx $scale -F+jTL -D${dx}/${dy} -O -K >> $ps
            awk '{print $2, $1}' $data_a | gmt psxy $rgn $scale_sub -W1.5p,blue,solid -O -K -t50 >> $ps
            awk '{print $3, $1}' $data_a | gmt psxy $rgn $scale_sub -W1.5p,red,solid -O -K -t50 >> $ps

            gmt psbasemap $rgn $scale_sub -Bpy${y_major_tick}f${y_minor_tick}+l"${y_label}" -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -X5c -Y0c -B${axis_b} -O -K >> $ps
            awk '{print $2, $1}' $data_b | gmt psxy $rgn $scale_sub -W1.5p,blue,solid -O -K -t50 >> $ps
            awk '{print $3, $1}' $data_b | gmt psxy $rgn $scale_sub -W1.5p,red,solid -O -K -t50 >> $ps
            dx=4.85
            dy=-0.2
            echo "0 1 ${depthtop_b}-${depthbottom_b} km" | gmt pstext $rgnx $scale -F+jTR -D${dx}/${dy} -O -K >> $ps

            idx_line=$(echo $idx_line | awk '{print $1+1}')

        done

    done

    gmt psbasemap $rgn $scale_sub -B0 -O >> $ps
    gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

plot_geodetic_vs_model_viscosity_kdes(){

    # set reasonable estimates of geodetic 'observed' viscosity (fitted to observed uplift rates - sig. trade-off with LAB depth)
    echo -n > junk
    echo "ASE B18 18.4 19.4 150 175 150 250 c -0c -12.5c NesW Nesw 0 3.5 B18" >> junk
    echo "S21AP S21 17.5 19.0 125 175 125 250 d 6c 0c Nesw NEsw 0 6.5 S21" >> junk
    echo "I11AP I11 19.3 20.0 125 175 125 250 e -16c -9c neSW neSw 0 2.5 I11" >> junk
    echo "Wo15AP Wo15 20.0 20.5 125 175 125 250 f 6c 0c neSw nESw 0 2.5 W15" >> junk

    # Set GMT plotting parameters   
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain FORMAT_GEO_OUT F FORMAT_GEO_MAP F MAP_ANNOT_OBLIQUE 15

    fold_plot_output_histogram=${fold_plot_output}
    mkdir -p ${fold_plot_output_histogram}
    plot_name="${fold_plot_output_histogram}/viscosity_kdes"
    ps="${plot_name}.ps"
    jpg="${plot_name}.jpg"

    rgnx="-R0/1/0/1"
    scale="-JX10c/8c"
    scale_sub="-JX5c/8c"

    rgn_map="-R-84/-65/-36/-70r"
    proj_map="-JS0/-90/10c"
    gmt psbasemap $rgn_map $proj_map -X30c -Y30c -B0 -K > $ps
    gmt makecpt -Chot -D -T18/23/0.1 -I -G0.08/1 > visc.cpt
    logetagrid=${fold_data_input_thermo}/${date}/Tp_${Tp}_sol50_${sol50}/viscosity_grids/distribution/summary/viscosity_depth_150_Tp_${Tp}_sol50_${sol50}_mean.grd
    gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A5000/0/4 -W1.2pt,black -O -K >> $ps
    awk '{if(NR>=4) print $2, $3, "100k"}' ${fold_data_input_locs}/sample_locs.txt |\
     gmt psxy $rgn_map $proj_map -SE- -W2p,blue,solid -O -K >> $ps
    for ((idx_loc=4; idx_loc<=6; idx_loc++)); do
        dx=$(awk '{if(NR=='${idx_loc}') print $5}' ${fold_data_input_locs}/sample_locs.txt)
        dy=$(awk '{if(NR=='${idx_loc}') print $6}' ${fold_data_input_locs}/sample_locs.txt)
        echo $idx_loc $dx $dy
        awk '{if(NR=='${idx_loc}') print $2, $3, $4}' ${fold_data_input_locs}/sample_locs.txt |\
        gmt pstext $rgn_map $proj_map -F+jBL -F+f18p,blue -D${dx}/${dy} -O -K >> $ps
    done
    dx=9.83
    dy=-9.75
    #gmt psbasemap $rgn_map $proj_map -Bxa5f2.5 -Bpya2f1 -BNESw -O -K >> $ps
    gmt psbasemap $rgn_map $proj_map -Bx15f1g15 -By5f1g5 -BNESw -O -K >> $ps
    echo "0 1 b" | gmt pstext $rgnx $proj_map -F+jBR -W1p,black,solid -Gwhite -TO -D${dx}/${dy} -O -K >> $ps
    gmt psscale -Dx-4.5c/-2.0c+w8c/0.25c+e+h -Cvisc.cpt -B1f0.5+l"log@-10@-@~\155@~@-@~\150@~@- (Pa s)" -Al -O -K >> $ps

    rgn_map="-R-117/-66/-88/-79.5r"
    proj_map="-JS0/-90/10c"
    gmt psbasemap $rgn_map $proj_map -X-11c -B0 -O -K >> $ps
    gmt grdimage $logetagrid -Cvisc.cpt $rgn_map $proj_map -O -K >> $ps
    gmt pscoast $rgn_map $proj_map -Dh -A5000/0/4 -W1.2pt,black -O -K >> $ps
    awk '{if(NR==2) print $2, $3, "100k"}' ${fold_data_input_locs}/sample_locs.txt |\
     gmt psxy $rgn_map $proj_map -SE- -W2p,blue,solid -O -K >> $ps
    dx=$(awk '{if(NR==2) print $5}' ${fold_data_input_locs}/sample_locs.txt)
    dy=$(awk '{if(NR==2) print $6}' ${fold_data_input_locs}/sample_locs.txt)
    awk '{if(NR==2) print $2, $3, $4}' ${fold_data_input_locs}/sample_locs.txt |\
     gmt pstext $rgn_map $proj_map -F+jBL -F+f18p,blue -D${dx}/${dy} -O -K >> $ps
    dx=0.15
    dy=-9.75
    gmt psbasemap $rgn_map $proj_map -Bx15f1g15 -By5f1g5 -BNeSW -O -K >> $ps
    echo "0 1 a" | gmt pstext $rgnx $proj_map -F+jBL -W1p,black,solid -Gwhite -TO -D${dx}/${dy} -O -K >> $ps

    idx_line=$(echo 1)
    idx_line_max=$(cat junk | wc -l)

    for ((idx_i=0; idx_i<=1; idx_i++)); do

        for ((idx_j=0; idx_j<=1; idx_j++)); do

            move_x=$(awk '{if(NR=='${idx_line}') print $10}' junk)
            move_y=$(awk '{if(NR=='${idx_line}') print $11}' junk)
            axis_a=$(awk '{if(NR=='${idx_line}') print $12}' junk)
            axis_b=$(awk '{if(NR=='${idx_line}') print $13}' junk)
            axis_x=$(echo $idx_i $idx_j | awk '{print "ns"}')
            axis_y=$(echo $idx_i $idx_j | awk '{if($2<1) print "We"; else print "wE"}')

            # set key variables incl. location, loading model, geodetic viscosity estimate, and file name for apparent viscosity data
            loc=$(awk '{if(NR=='${idx_line}') print $1}' junk)
            loading_model=$(awk '{if(NR=='${idx_line}') print $2}' junk)
            loading_model_label=$(awk '{if(NR=='${idx_line}') print $16}' junk)
            echo $loc $loading_model $idx_i $idx_j
            eta_min=$(awk '{if(NR=='${idx_line}') print $3}' junk)
            eta_max=$(awk '{if(NR=='${idx_line}') print $4}' junk)

            depthtop_a=$(awk '{if(NR=='${idx_line}') print $5}' junk)
            depthbottom_a=$(awk '{if(NR=='${idx_line}') print $6}' junk)
            data_a=${fold_data_output}/apparent_viscosity/distribution/summary/apparent_viscosity_location_${loc}_loading_${loading_model}_depthtop_${depthtop_a}_depthbottom_${depthbottom_a}_kde.txt
            depthtop_b=$(awk '{if(NR=='${idx_line}') print $7}' junk)
            depthbottom_b=$(awk '{if(NR=='${idx_line}') print $8}' junk)
            data_b=${fold_data_output}/apparent_viscosity/distribution/summary/apparent_viscosity_location_${loc}_loading_${loading_model}_depthtop_${depthtop_b}_depthbottom_${depthbottom_b}_kde.txt
            label=$(awk '{if(NR=='${idx_line}') print $9}' junk)

            # set histogram plot region and label
            y_min=17.5
            y_max=21.0
            y_major_tick=1
            y_minor_tick=0.1
            bin_width=0.05
            x_min=$(awk '{if(NR=='${idx_line}') print $14}' junk)
            x_max=$(awk '{if(NR=='${idx_line}') print $15}' junk)
            x_major_tick=1
            x_minor_tick=0.5
            x_label="f(log@-10@-@~\150@~)"
            y_label="log@-10@-@~\150@~ (Pa s)"
            rgn="-R${x_min}/${x_max}/${y_min}/${y_max}"
            
            # plot histogram
            gmt psbasemap $rgn $scale_sub -Bpy${y_major_tick}f${y_minor_tick}+l"${y_label}" -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -X$move_x -Y$move_y -B${axis_a} -O -K >> $ps
            echo 0 ${eta_min} > junk2
            echo 0 ${eta_max} >> junk2
            echo ${x_max} ${eta_max} >> junk2
            echo ${x_max} ${eta_min} >> junk2
            head -n 1 junk2 >> junk2
            awk '{print $1, $2}' junk2 | gmt psxy $rgn $scale -Ggray -t50 -O -K >> $ps
            dx=9.85
            dy=-7.8
            echo "0 1 ${loading_model_label}" | gmt pstext $rgnx $scale -F+jBR -D${dx}/${dy} -O -K >> $ps
            dx=0.15
            dy=-7.8
            echo "0 1 ${label}" | gmt pstext $rgnx $scale -F+jBL -W1p,black,solid -Gwhite -TO -D${dx}/${dy} -O -K >> $ps
            dx=0.15
            dy=-0.2
            echo "0 1 ${depthtop_a}-${depthbottom_a} km" | gmt pstext $rgnx $scale -F+jTL -D${dx}/${dy} -O -K >> $ps
            awk '{print $2, $1}' $data_a | gmt psxy $rgn $scale_sub -W1.5p,blue,solid -O -K -t50 >> $ps
            awk '{print $3, $1}' $data_a | gmt psxy $rgn $scale_sub -W1.5p,red,solid -O -K -t50 >> $ps

            gmt psbasemap $rgn $scale_sub -Bpy${y_major_tick}f${y_minor_tick}+l"${y_label}" -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -X5c -Y0c -B${axis_b} -O -K >> $ps
            awk '{print $2, $1}' $data_b | gmt psxy $rgn $scale_sub -W1.5p,blue,solid -O -K -t50 >> $ps
            awk '{print $3, $1}' $data_b | gmt psxy $rgn $scale_sub -W1.5p,red,solid -O -K -t50 >> $ps
            dx=4.85
            dy=-0.2
            echo "0 1 ${depthtop_b}-${depthbottom_b} km" | gmt pstext $rgnx $scale -F+jTR -D${dx}/${dy} -O -K >> $ps

            idx_line=$(echo $idx_line | awk '{print $1+1}')

        done

    done

    gmt psbasemap $rgn $scale_sub -B0 -O >> $ps
    gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

plot_geodetic_vs_model_viscosity_distributions(){

    # set reasonable estimates of geodetic 'observed' viscosity (fitted to observed uplift rates - sig. trade-off with LAB depth)
    echo -n > junk
    echo "AP S20 17.7 19.0" >> junk
    echo "AP Wo15 20.0 20.5" >> junk
    echo "AP I11 19.3 20.0" >> junk
    echo "ASE B18 18.6 19.4" >> junk
    
    # Set GMT plotting parameters   
    gmt gmtset FORMAT_FLOAT_OUT %.8f FONT_LABEL 16p FONT_ANNOT_PRIMARY 14p PS_MEDIA a0 MAP_FRAME_TYPE plain

    fold_plot_output_histogram=${fold_plot_output}
    mkdir -p ${fold_plot_output_histogram}
    plot_name="${fold_plot_output_histogram}/viscosity_distributions"
    ps="${plot_name}.ps"
    jpg="${plot_name}.jpg"

    rgn="-R0/1/0/1"
    scale="-JX10c/10c"

    gmt psbasemap $rgn $scale -X30c -Y30c -B0 -K > $ps

    idx_line=1
    idx_line_max=$(cat junk | wc -l)

    # set key variables incl. location, loading model, geodetic viscosity estimate, and file name for apparent viscosity data
    loc=$(awk '{if(NR=='${idx_line}') print $1}' junk)
    loading_model=$(awk '{if(NR=='${idx_line}') print $2}' junk)
    eta_min=$(awk '{if(NR=='${idx_line}') print $3}' junk)
    eta_max=$(awk '{if(NR=='${idx_line}') print $4}' junk)
    data=${fold_data_output}/apparent_viscosity/distribution/individual/apparent_viscosity_location_${loc}_loading_${loading_model}_all.txt
    
    # set histogram plot region and label
    x_min=17
    x_max=22
    x_major_tick=0.5
    x_minor_tick=0.1
    bin_width=0.05
    y_max=25
    x_label="log@-10@-@~\150@~ (location ${loc}, loading ${loading_model})"
    rgn="-R${x_min}/${x_max}/0/${y_max}"
    
    # plot histogram
    gmt psbasemap $rgn $scale -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -Bpy10f5+l"Frequency (%)" -BWSne -O -K >> $ps
    awk '{print $1}' $data | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gblue -O -K -t50 >> $ps
    awk '{print $2}' $data | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gred -O -K -t50 >> $ps
    echo ${eta_min} 0 > junk2
    echo ${eta_max} 0 >> junk2
    echo ${eta_max} ${y_max} >> junk2
    echo ${eta_min} ${y_max} >> junk2
    head -n 1 junk2 >> junk2
    awk '{print $1, $2}' junk2 | gmt psxy $rgn $scale -Ggray -t50 -O -K >> $ps

    for ((idx_line=2; idx_line<=idx_line_max; idx_line++)); do

        # set key variables incl. location, loading model, geodetic viscosity estimate, and file name for apparent viscosity data
        loc=$(awk '{if(NR=='${idx_line}') print $1}' junk)
        loading_model=$(awk '{if(NR=='${idx_line}') print $2}' junk)
        eta_min=$(awk '{if(NR=='${idx_line}') print $3}' junk)
        eta_max=$(awk '{if(NR=='${idx_line}') print $4}' junk)
        data=${fold_data_output}/apparent_viscosity/distribution/individual/apparent_viscosity_location_${loc}_loading_${loading_model}_all.txt
        
        # set histogram plot region and label
        x_min=17
        x_max=22
        x_major_tick=0.5
        x_minor_tick=0.1
        bin_width=0.05
        y_max=25
        x_label="log@-10@-@~\150@~ (location ${loc}, loading ${loading_model})"
        rgn="-R${x_min}/${x_max}/0/${y_max}"
        move="-X11c"
        
        # plot histogram
        gmt psbasemap $rgn $scale $move -Bpx${x_major_tick}f${x_minor_tick}+l"${x_label}" -Bpy10f5+l"Frequency (%)" -BwSne -O -K >> $ps
        awk '{print $1}' $data | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gblue -O -K -t50 >> $ps
        awk '{print $2}' $data | gmt pshistogram -Z1 $rgn $scale -W${bin_width} -Gred -O -K -t50 >> $ps
        rm -f junk2
        echo ${eta_min} 0 > junk2
        echo ${eta_max} 0 >> junk2
        echo ${eta_max} ${y_max} >> junk2
        echo ${eta_min} ${y_max} >> junk2
        head -n 1 junk2 >> junk2
        awk '{print $1, $2}' junk2 | gmt psxy $rgn $scale -Ggray -t50 -O -K >> $ps

    done

    gmt psbasemap $rgn $scale -B0 -O >> $ps
    gmt psconvert $ps -Tj -E600 -A0.1c -P -Z

}

fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')
input_locs=$(awk '$1 ~ /^input_locations/' config.ini | awk '{print $3}')
fold_data_input_locs=${fold_data_input}/${input_locs}
# since this is a data preparation stage, data output will be routed towards an input dir
Tp=$(awk '$1 ~ /^potential_temperature/' config.ini | awk '{print $3}')
sol50=$(awk '$1 ~ /^solidus_50km/' config.ini | awk '{print $3}')
date=$(awk '$1 ~ /^date/' config.ini | awk '{print $3}')
fold_data_output=${fold_data_output}/${date}/Tp_${Tp}_sol50_${sol50}
fold_plot_output=${fold_plot_output}/${date}/Tp_${Tp}_sol50_${sol50}
fold_data_input_viscosity=$(awk '$1 ~ /^input_viscosity/' config.ini | awk '{print $3}')
input_viscosity=${fold_data_input}/${fold_data_input_viscosity}/${date}/Tp_${Tp}_sol50_${sol50}
fold_data_input_thermo=${fold_data_input}/$(awk '$1 ~ /^input_thermodynamic/' config.ini | awk '{print $3}')

mkdir -p ${fold_plot_output}

plot_geodetic_vs_model_viscosity_kdes
#plot_geodetic_vs_model_viscosity_kdes_ellipses

rm -f junk*