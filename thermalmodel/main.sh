#!/bin/bash

get_melting_profile(){

    echo "Determining melting profile for potential temperature = $Tp degC"

    echo -n > ${meltfold}/melt_output_${type}_${Tp}.dat
    ${meltfold}/3phmlt.o -S -x0 -p1 -h0 -T${Tp} >> ${meltfold}/melt_output_${type}_${Tp}.dat 2>&1
    ####################################################################
    ## Convert pressure into depth
    ####################################################################
    awk 'NR>6{print $1, $2, $6}' ${meltfold}/melt_output_${type}_${Tp}.dat | sort -nk1 \
        > ${meltfold}/melting_region_solidus_${type}_${Tp}.PTF

    minp=$(awk 'NR==1{print $1}' ${meltfold}/melting_region_solidus_${type}_${Tp}.PTF)
    crust=$(gmt gmtmath -Q $minp 10 6 POW MUL ${rhoc} 9.81 MUL DIV =)
    echo "crustal thickness is $crust for Tp of $Tp"
    minp_k=$(gmt gmtmath -Q ${rhoc} 9.81 MUL 7000 MUL 1e9 DIV =)
    awk '{print (('$minp_k')*(10^6)/('${rhoc}'*9.81))+(($1-'$minp')*30), $2, $3}' ${meltfold}/melting_region_solidus_${type}_${Tp}.PTF \
        > tmp
    T1=$(awk 'NR==1{print $2}' tmp)
    T2=$(awk 'NR==2{print $2}' tmp)
    Z1=$(awk 'NR==1{print $1}' tmp)
    Z2=$(awk 'NR==2{print $1}' tmp)
    F1=$(awk 'NR==1{print $3}' tmp)
    grad=$(echo $Z1 $Z2 $T1 $T2 | awk '{print ($4-$3)/($2-$1)}')
    echo "0 0 $F1" > header_crust
    cat header_crust tmp | awk '{print $2, $3, $1}' | uniq -f 2 | awk '{print $3, $1, $2}' \
        | gmt sample1d -Fl -I1 > ${meltfold}/melting_region_solidus_crust_${type}_${Tp}.zTF
    awk '{print $2}' ${meltfold}/melting_region_solidus_crust_${type}_${Tp}.zTF > ${codefold}/initial_geotherms/mac${Tp}.dat
    rm -f tmp header*

}


plate_model() {

    ####################################################################
    # Set up plate model
    ####################################################################
    echo "Running plate model for temperature = $Tp degC; plate thickness = $zp km; zero-age ridge depth = $rd m"

    cd ${codefold}
    if [ $Tp -eq $Tpmin -a $zp -eq $zpmin -a $rd -eq $rdmin ]; then
        make clean
        make
        ./cool ${Tp} ${zp} ${rd}
    else
        rm -f ${outfold}/Tgrid-${Tp}-${zp}-${rd}.dat
        ./cool ${Tp} ${zp} ${rd}
    fi
    if [ $temp_flag -eq 0 ]; then
    	rm -f ${outfold}/Tgrid-${Tp}-${zp}-${rd}.dat
    elif [ $temp_flag -eq 1 ]; then
        echo "Keeping temperature file"
    else
        echo "Unexpected value for temperature file flag"
    fi
    cd ..

}

assess_joint_misfit(){

    ###################################################################
    # Subsidence
    ###################################################################
    echo "Determining subsidence misfit for $Tp degC"
    if [ $Tp -eq $Tpmin ]; then
        echo -n > ${outfold}/RMS_misfit_depth_${region}.txt
    fi
    awk '{if ($4>=5 && $4<=200){print $4}}' ${file_ref_subs} | sort -nk1 > xpoints
    awk '{if ($4>=5 && $4<=200){print $4, -($3*1000), '${sigma_subs}'}}' ${file_ref_subs} | sort -k1,1n > input
    for zp in $(seq $zpmin $zpinc $zpmax); do
        for rd in $(seq $rdmin $rdinc $rdmax); do
            echo "Plate thickness is $zp km, ridge depth is $rd"
            file="${outfold}/depth-${Tp}-${zp}-${rd}.dat"
            gmt sample1d ${file} -Nxpoints -T0 -Fa | awk '{print -$2}' > resampled
            paste input resampled | awk '{sum=sum+((($2-$4)/$3)^2)}END{print '$Tp', \
                '$zp', '$rd', sqrt(sum/NR)}' >> ${outfold}/RMS_misfit_depth_${region}.txt
            tail -n1 ${outfold}/RMS_misfit_depth_${region}.txt
        done
    done

    ###################################################################
    # Heat Flow
    ###################################################################
    echo "Determining heatflow misfit for $Tp degC"
    if [ $Tp -eq $Tpmin ]; then
        echo -n > ${outfold}/RMS_misfit_heatflow_${region}.txt
    fi
    awk '{if (NR>1 && $1>=5 && $1<=168){print $1}}' ${file_ref_hf} | sort -nk1 > xpoints
    awk '{if (NR>1 && $1>=5 && $1<=168){print $1, $4, $6-$5}}' ${file_ref_hf} | sort -nk1 > input
    for zp in $(seq $zpmin $zpinc $zpmax); do
        for rd in $(seq $rdmin $rdinc $rdmax); do
            echo "Plate thickness is $zp km, ridge depth is $rd"
            file="${outfold}/hf-${Tp}-${zp}-${rd}.dat"
            awk '{if ($1>=5){print $0}}' ${file} | gmt sample1d -Nxpoints -T0 -Fa | awk '{print $2}' > resampled
            paste input resampled | awk '{sum=sum+(((1.349*($2-$4))/$3)^2)}END{print '$Tp','$zp', '$rd', sqrt(sum/NR)}' >> \
                ${outfold}/RMS_misfit_heatflow_${region}.txt
        done
    done

    ###################################################################
    # Joint Misfit
    ###################################################################
    echo "Determining total misfit for $Tp degC"
    if [ $Tp -eq $Tpmin ]; then
        echo -n > ${outfold}/RMS_misfit_total_${region}.txt
    fi
    for zp in $(seq $zpmin $zpinc $zpmax); do
        for rd in $(seq $rdmin $rdinc $rdmax); do
            Ms=$(awk '{if ($1=='$Tp' && $2=='$zp' && $3=='$rd'){print $4}}' ${outfold}/RMS_misfit_depth_${region}.txt)
            Mq=$(awk '{if ($1=='$Tp' && $2=='$zp' && $3=='$rd'){print $4}}' ${outfold}/RMS_misfit_heatflow_${region}.txt)
            Mt=$(echo $Ms $Mq | awk '{print sqrt(((($1*'${weight_subs}')^2)+(($2*'${weight_heatflow}')^2))/2)}')
            echo $Tp $zp $rd $Mt >> ${outfold}/RMS_misfit_total_${region}.txt
        done
    done

}

assess_joint_misfit_singleT(){

    ###################################################################
    # Subsidence
    ###################################################################
    echo "Determining subsidence misfit for $Tp degC"
    echo -n > ${outfold}/RMS_misfit_depth_${Tp}_${region}.txt
    file_ref="${subsfold}/${region}_basementdepth_1d.dat"
    awk '{if ($4>=5 && $4<=200){print $4}}' ${file_ref_subs} | sort -nk1 > xpoints
    awk '{if ($4>=5 && $4<=200){print $4, -($3*1000), '${sigma_subs}'}}' ${file_ref_subs} | sort -k1,1n > input
    for zp in $(seq $zpmin $zpinc $zpmax); do
        for rd in $(seq $rdmin $rdinc $rdmax); do
            echo "Plate thickness is $zp km, ridge depth is $rd"
            file="${outfold}/depth-${Tp}-${zp}-${rd}.dat"
            gmt sample1d ${file} -Nxpoints -T0 -Fa | awk '{print -$2}' > resampled
            paste input resampled | awk '{sum=sum+((($2-$4)/$3)^2)}END{print '$Tp', \
                '$zp', '$rd', sqrt(sum/NR)}' >> ${outfold}/RMS_misfit_depth_${Tp}_${region}.txt
        done
    done

    ###################################################################
    # Heat Flow
    ###################################################################
    echo "Determining heatflow misfit for $Tp degC"
    echo -n > ${outfold}/RMS_misfit_heatflow_${Tp}_${region}.txt
    awk '{if (NR>1 && $1>=5 && $1<=168){print $1}}' ${file_ref_hf} | sort -nk1 > xpoints
    awk '{if (NR>1 && $1>=5 && $1<=168){print $1, $4, $6-$5}}' ${file_ref_hf} | sort -nk1 > input
    for zp in $(seq $zpmin $zpinc $zpmax); do
        for rd in $(seq $rdmin $rdinc $rdmax); do
            echo "Plate thickness is $zp km, ridge depth is $rd"
            file="${outfold}/hf-${Tp}-${zp}-${rd}.dat"
            awk '{if ($1>=5){print $0}}' ${file} | gmt sample1d -Nxpoints -T0 -Fa | awk '{print $2}' > resampled
            paste input resampled | awk '{sum=sum+(((1.349*($2-$4))/$3)^2)}END{print '$Tp','$zp', '$rd', sqrt(sum/NR)}' >> \
                ${outfold}/RMS_misfit_heatflow_${Tp}_${region}.txt
        done
    done

    ###################################################################
    # Joint Misfit
    ###################################################################
    echo "Determining total misfit for $Tp degC"
    echo -n > ${outfold}/RMS_misfit_total_${Tp}_${region}.txt
    for zp in $(seq $zpmin $zpinc $zpmax); do
        for rd in $(seq $rdmin $rdinc $rdmax); do
            Ms=$(awk '{if ($1=='$Tp' && $2=='$zp' && $3=='$rd'){print $4}}' ${outfold}/RMS_misfit_depth_${Tp}_${region}.txt)
            Mq=$(awk '{if ($1=='$Tp' && $2=='$zp' && $3=='$rd'){print $4}}' ${outfold}/RMS_misfit_heatflow_${Tp}_${region}.txt)
            Mt=$(echo $Ms $Mq | awk '{print sqrt(((($1*'${weight_subs}')^2)+(($2*'${weight_heatflow}')^2))/2)}')
            echo $Tp $zp $rd $Mt >> ${outfold}/RMS_misfit_total_${Tp}_${region}.txt
        done
    done

}

plot_all_misfit(){

    ps="${fold_plot_output}/all_misfit_${region}.ps"
    jpg="${fold_plot_output}/all_misfit_${region}.jpg"


    Tpbest=$(awk 'NR==1 || $4 < min {min=$4; Tp=$1}END{print Tp}' ${outfold}/RMS_misfit_depth_${region}.txt)
    zpbest=$(awk 'NR==1 || $4 < min {min=$4; zp=$2}END{print zp}' ${outfold}/RMS_misfit_depth_${region}.txt)
    rdbest=$(awk 'NR==1 || $4 < min {min=$4; rd=$3}END{print rd}' ${outfold}/RMS_misfit_depth_${region}.txt)
    echo "Best params for depth only: rd=$rdbest, Tp=$Tpbest, zp=$zpbest"

    inc_Tprd="-I1/4"
    inc_rdzp="-I4/0.4"
    inc_Tpzp="-I1/0.4"

    awk '{if ($2=='$zpbest'){print $1, $3, $4}}' ${outfold}/RMS_misfit_depth_${region}.txt \
        | gmt surface $inc_Tprd -G${outfold}/Tprd_misfit_depth_${region}.grd -R$Tpmin/$Tpmax/$rdmin/$rdmax
    awk '{if ($1=='$Tpbest'){print $3, $2, $4}}' ${outfold}/RMS_misfit_depth_${region}.txt \
        | gmt surface $inc_rdzp -G${outfold}/rdzp_misfit_depth_${region}.grd -R$rdmin/$rdmax/$zpmin/$zpmax

    # Plot RMS plate thickness vs temperature:
    awk '{if ($3=='$rdbest'){print $1, $2, $4}}' ${outfold}/RMS_misfit_depth_${region}.txt | gmt blockmean -I1/0.4 \
        -R$Tpmin/$Tpmax/$zpmin/$zpmax | gmt surface -I1/0.4 -G${outfold}/Tpzp_misfit_depth_${region}.grd \
            -R$Tpmin/$Tpmax/$zpmin/$zpmax

    gmt makecpt -Cgray -T0.5/1.5/0.01 -I -D > misfit.cpt
    gmt psbasemap -JX21/11.5c -R0/1/0/1 -B+t"${label_arr[$m]}" -K -Y20c > $ps
    gmt grdimage ${outfold}/Tpzp_misfit_depth_${region}.grd -JX10c -R$Tpmin/$Tpmax/$zpmin/$zpmax -Cmisfit.cpt -K -O \
        -Bx100f${rdinc}+l"Temperature (@~\260@~C)" -By10f${zpinc}+l"Plate Thickness (km)" -BWsNe >> $ps
    gmt grdcontour ${outfold}/Tpzp_misfit_depth_${region}.grd -R -J -C0.1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/Tpzp_misfit_depth_${region}.grd -R -J -A0.2 -Gn1 -K -O -Wthin >> $ps
    echo "1107.5 54 @~\143@~@+2@+@-s@-" | gmt pstext -R -J -F+f12p+jLB -TO -Gwhite -W1p -O -K >> $ps
    awk '{if ($3=='$rdbest'){print $1, $2, $4}}' ${outfold}/RMS_misfit_depth_${region}.txt | \
        awk 'NR==1 || $3<min {min=$3;Tp=$1;rd=$2}END{print Tp, rd}' | gmt psxy -R -J -Sc0.25c -Gred -Wthin -O -K >> $ps

    Tpbest=$(awk 'NR==1 || $4 < min {min=$4; Tp=$1}END{print Tp}' ${outfold}/RMS_misfit_heatflow_${region}.txt)
    zpbest=$(awk 'NR==1 || $4 < min {min=$4; zp=$2}END{print zp}' ${outfold}/RMS_misfit_heatflow_${region}.txt)
    rdbest=$(awk 'NR==1 || $4 < min {min=$4; rd=$3}END{print rd}' ${outfold}/RMS_misfit_heatflow_${region}.txt)
    echo "Best params for heatflow only: rd=$rdbest, Tp=$Tpbest, zp=$zpbest"

    awk '{if ($2=='$zpbest'){print $1, $3, $4}}' ${outfold}/RMS_misfit_heatflow_${region}.txt \
        | gmt surface $inc_Tprd -G${outfold}/Tprd_misfit_heatflow_${region}.grd -R$Tpmin/$Tpmax/$rdmin/$rdmax
    awk '{if ($1=='$Tpbest'){print $3, $2, $4}}' ${outfold}/RMS_misfit_heatflow_${region}.txt \
        | gmt surface $inc_rdzp -G${outfold}/rdzp_misfit_heatflow_${region}.grd -R$rdmin/$rdmax/$zpmin/$zpmax

    # Plot RMS plate thickness vs temperature:
    awk '{if ($3=='$rdbest'){print $1, $2, $4}}' ${outfold}/RMS_misfit_heatflow_${region}.txt | gmt blockmean -I1/0.4 \
        -R$Tpmin/$Tpmax/$zpmin/$zpmax | gmt surface -I1/0.4 -G${outfold}/Tpzp_misfit_heatflow_${region}.grd \
            -R$Tpmin/$Tpmax/$zpmin/$zpmax

    gmt grdimage ${outfold}/Tpzp_misfit_heatflow_${region}.grd -JX10c -R$Tpmin/$Tpmax/$zpmin/$zpmax -Cmisfit.cpt -K -O \
        -Bx100f${rdinc}+l"Temperature (@~\260@~C)" -By10f${zpinc}+l"Plate Thickness (km)" -BwsNE -X11c >> $ps
    gmt grdcontour ${outfold}/Tpzp_misfit_heatflow_${region}.grd -R -J -C0.1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/Tpzp_misfit_heatflow_${region}.grd -R -J -A0.2 -Gn1 -K -O -Wthin >> $ps
    echo "1107.5 54 @~\143@~@+2@+@-hf@-" | gmt pstext -R -J -F+f12p+jLB -TO -Gwhite -W1p -O -K >> $ps
    awk '{if ($3=='$rdbest'){print $1, $2, $4}}' ${outfold}/RMS_misfit_heatflow_${region}.txt | \
        awk 'NR==1 || $3<min {min=$3;Tp=$1;rd=$2}END{print Tp, rd}' | gmt psxy -R -J -Sc0.25c -Gred -Wthin -O -K >> $ps

    Tpbest=$(awk 'NR==1 || $4 < min {min=$4; Tp=$1}END{print Tp}' ${outfold}/RMS_misfit_total_${region}.txt)
    zpbest=$(awk 'NR==1 || $4 < min {min=$4; zp=$2}END{print zp}' ${outfold}/RMS_misfit_total_${region}.txt)
    rdbest=$(awk 'NR==1 || $4 < min {min=$4; rd=$3}END{print rd}' ${outfold}/RMS_misfit_total_${region}.txt)
    echo "Best params for total only: rd=$rdbest, Tp=$Tpbest, zp=$zpbest"

    awk '{if ($2=='$zpbest'){print $1, $3, $4}}' ${outfold}/RMS_misfit_total_${region}.txt \
        | gmt surface $inc_Tprd -G${outfold}/Tprd_misfit_total_${region}.grd -R$Tpmin/$Tpmax/$rdmin/$rdmax
    awk '{if ($1=='$Tpbest'){print $3, $2, $4}}' ${outfold}/RMS_misfit_total_${region}.txt \
        | gmt surface $inc_rdzp -G${outfold}/rdzp_misfit_total_${region}.grd -R$rdmin/$rdmax/$zpmin/$zpmax

    # Plot RMS plate thickness vs temperature:
    awk '{if ($3=='$rdbest'){print $1, $2, $4}}' ${outfold}/RMS_misfit_total_${region}.txt | gmt blockmean -I1/0.4 \
        -R$Tpmin/$Tpmax/$zpmin/$zpmax | gmt surface -I1/0.4 -G${outfold}/Tpzp_misfit_total_${region}.grd \
            -R$Tpmin/$Tpmax/$zpmin/$zpmax

    gmt grdimage ${outfold}/Tpzp_misfit_total_${region}.grd -JX10c -R$Tpmin/$Tpmax/$zpmin/$zpmax -Cmisfit.cpt -K -O \
        -Bx100f${rdinc}+l"Temperature (@~\260@~C)" -By10f${zpinc}+l"Plate Thickness (km)" -BWSnE -X-5.5c -Y-11c >> $ps
    gmt grdcontour ${outfold}/Tpzp_misfit_total_${region}.grd -R -J -C0.1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/Tpzp_misfit_total_${region}.grd -R -J -A0.2 -Gn1 -K -O -Wthin >> $ps
    echo "1107.5 54 @~\143@~@+2@+@-t@-" | gmt pstext -R -J -F+f12p+jLB -TO -Gwhite -W1p -O -K >> $ps
    awk '{if ($3=='$rdbest'){print $1, $2, $4}}' ${outfold}/RMS_misfit_total_${region}.txt | \
        awk 'NR==1 || $3<min {min=$3;Tp=$1;rd=$2}END{print Tp, rd}' | gmt psxy -R -J -Sc0.25c -Gred -Wthin -O >> $ps
    gmt psconvert -P $ps -Tj -E300 -A0.5c/0.5c -Z


}

plot_singleT_misfit(){

    ps="${fold_plot_output}/all_misfit_${Tp}_${region}.ps"
    jpg="${fold_plot_output}/all_misfit_${Tp}_${region}.jpg"

    zpbest=$(awk 'NR==1 || $4 < min {min=$4; zp=$2}END{print zp}' ${outfold}/RMS_misfit_depth_${Tp}_${region}.txt)
    rdbest=$(awk 'NR==1 || $4 < min {min=$4; rd=$3}END{print rd}' ${outfold}/RMS_misfit_depth_${Tp}_${region}.txt)
    echo "Best params for old depth only: rd=$rdbest, Tp=$Tp, zp=$zpbest"

    # Plot RMS plate thickness vs temperature:
    awk '{print $3, $2, $4}' ${outfold}/RMS_misfit_depth_${Tp}_${region}.txt | gmt blockmean -I4/0.4 \
        -R$rdmin/$rdmax/$zpmin/$zpmax | gmt surface -I4/0.4 -G${outfold}/rdzp_misfit_depth_${Tp}_${region}.grd \
            -R$rdmin/$rdmax/$zpmin/$zpmax

    gmt makecpt -Cgray -T0.5/1.5/0.01 -I -D > misfit.cpt

    ref=1.05
    cont=$(gmt grdinfo ${outfold}/rdzp_misfit_depth_${Tp}_${region}.grd -C | awk '{print $6*'$ref'}')
    gmt grdcontour ${outfold}/rdzp_misfit_depth_${Tp}_${region}.grd -R -J -C+$cont -D0.05_contour_rdzp_${Tp}_${region}.txt
    gmt grdimage ${outfold}/rdzp_misfit_depth_${Tp}_${region}.grd -JX10c -R$rdmin/$rdmax/$zpmin/$zpmax -Cmisfit.cpt -K \
        -Bx200f${rdinc}+l"Ridge Depth (m)" -By10f${zpinc}+l"Plate Thickness (km)" -BWNse -Y25c > $ps
    gmt grdcontour ${outfold}/rdzp_misfit_depth_${Tp}_${region}.grd -R -J -C0.1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_depth_${Tp}_${region}.grd -R -J -A0.2 -Gn1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_depth_${Tp}_${region}.grd -R -J -C+$cont -K -O -Wthin,black,'-' >> $ps
    echo "2020 54 @~\143@~@+2@+@-t@-" | gmt pstext -R -J -F+f12p+jLB -TO -Gwhite -W1p -O -K >> $ps
    awk '{print $3, $2, $4}' ${outfold}/RMS_misfit_depth_${Tp}_${region}.txt | \
        awk 'NR==1 || $3<min {min=$3;rd=$1;zp=$2}END{print rd, zp}' | gmt psxy -R -J -Sc0.25c -Gred -Wthin -O -K >> $ps

    zpbest=$(awk 'NR==1 || $3 < min {min=$3; zp=$2}END{print zp}' ${outfold}/RMS_misfit_heatflow_${Tp}_${region}.txt)
    echo "Best params for heatflow only: rd=$rdbest, Tp=$Tp, zp=$zpbest"

    # Plot RMS plate thickness vs temperature:
    awk '{print $3, $2, $4}' ${outfold}/RMS_misfit_heatflow_${Tp}_${region}.txt | gmt blockmean -I4/0.4 \
        -R$rdmin/$rdmax/$zpmin/$zpmax | gmt surface -I4/0.4 -G${outfold}/rdzp_misfit_heatflow_${Tp}_${region}.grd \
            -R$rdmin/$rdmax/$zpmin/$zpmax
    cont=$(gmt grdinfo ${outfold}/rdzp_misfit_heatflow_${Tp}_${region}.grd -C | awk '{print $6*'$ref'}')
    gmt grdcontour ${outfold}/rdzp_misfit_heatflow_${Tp}_${region}.grd -R -J -C+$cont -D0.05_contour_rdzp_${Tp}_${region}.txt
    gmt grdimage ${outfold}/rdzp_misfit_heatflow_${Tp}_${region}.grd -JX10c -R$rdmin/$rdmax/$zpmin/$zpmax -Cmisfit.cpt -K -O \
        -Bx200f${rdinc}+l"Ridge Depth (m)" -By10f${zpinc}+l"Plate Thickness (km)" -BwsNE -X11c >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_heatflow_${Tp}_${region}.grd -R -J -C0.1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_heatflow_${Tp}_${region}.grd -R -J -A0.2 -Gn1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_heatflow_${Tp}_${region}.grd -R -J -C+$cont -K -O -Wthin,black,'-' >> $ps
    echo "2020 54 @~\143@~@+2@+@-hf@-" | gmt pstext -R -J -F+f12p+jLB -TO -Gwhite -W1p -O -K >> $ps
    awk '{print $3, $2, $4}' ${outfold}/RMS_misfit_heatflow_${Tp}_${region}.txt | \
        awk 'NR==1 || $3<min {min=$3;rd=$1;zp=$2}END{print rd, zp}' | gmt psxy -R -J -Sc0.25c -Gred -Wthin -O -K >> $ps

    zpbest=$(awk 'NR==1 || $4 < min {min=$4; zp=$2}END{print zp}' ${outfold}/RMS_misfit_total_${Tp}_${region}.txt)
    rdbest=$(awk 'NR==1 || $4 < min {min=$4; rd=$3}END{print rd}' ${outfold}/RMS_misfit_total_${Tp}_${region}.txt)
    echo "Best params for total only: rd=$rdbest, Tp=$Tp, zp=$zpbest"
    # Plot RMS plate thickness vs temperature:
    awk '{print $3, $2, $4}' ${outfold}/RMS_misfit_total_${Tp}_${region}.txt | gmt blockmean -I4/0.4 \
        -R$rdmin/$rdmax/$zpmin/$zpmax | gmt surface -I4/0.4 -G${outfold}/rdzp_misfit_total_${Tp}_${region}.grd \
            -R$rdmin/$rdmax/$zpmin/$zpmax
    cont=$(gmt grdinfo ${outfold}/rdzp_misfit_total_${Tp}_${region}.grd -C | awk '{print $6*'$ref'}')
    gmt grdcontour ${outfold}/rdzp_misfit_total_${Tp}_${region}.grd -R -J -C+$cont -D0.05_contour_rdzp_${Tp}_${region}.txt
    gmt grdimage ${outfold}/rdzp_misfit_total_${Tp}_${region}.grd -JX10c -R$rdmin/$rdmax/$zpmin/$zpmax -Cmisfit.cpt -K -O \
        -Bx200f${rdinc}+l"Ridge Depth (m)" -By10f${zpinc}+l"Plate Thickness (km)" -BWSnE -X-5.5c -Y-11c >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_total_${Tp}_${region}.grd -R -J -C0.1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_total_${Tp}_${region}.grd -R -J -A0.2 -Gn1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_total_${Tp}_${region}.grd -R -J -C+$cont -K -O -Wthin,black,'-' >> $ps
    echo "2020 54 @~\143@~@+2@+@-t@-" | gmt pstext -R -J -F+f12p+jLB -TO -Gwhite -W1p -O -K >> $ps
    awk '{print $3, $2, $4}' ${outfold}/RMS_misfit_total_${Tp}_${region}.txt | \
        awk 'NR==1 || $3<min {min=$3;rd=$1;zp=$2}END{print rd, zp}' | gmt psxy -R -J -Sc0.25c -Gred -Wthin -O >> $ps
    gmt psconvert -P $ps -Tj -E300 -A0.5c/0.5c -Z



}


plot_joint_misfit(){

    ps="${fold_plot_output}/all_joint_misfit_${region}.ps"
    jpg="${fold_plot_output}/all_joint_misfit_${region}.jpg"

    Tpbest=$(awk 'NR==1 || $4 < min {min=$4; Tp=$1}END{print Tp}' ${outfold}/RMS_misfit_total_${region}.txt)
    zpbest=$(awk 'NR==1 || $4 < min {min=$4; zp=$2}END{print zp}' ${outfold}/RMS_misfit_total_${region}.txt)
    rdbest=$(awk 'NR==1 || $4 < min {min=$4; rd=$3}END{print rd}' ${outfold}/RMS_misfit_total_${region}.txt)

    inc_Tprd="-I1/4"
    inc_rdzp="-I4/0.4"
    inc_Tpzp="-I1/0.4"

    awk '{if ($2=='$zpbest'){print $1, $3, $4}}' ${outfold}/RMS_misfit_total_${region}.txt \
        | gmt surface $inc_Tprd -G${outfold}/Tprd_misfit_total_${region}.grd \
            -R$Tpmin/$Tpmax/$rdmin/$rdmax
    awk '{if ($1=='$Tpbest'){print $3, $2, $4}}' ${outfold}/RMS_misfit_total_${region}.txt \
        | gmt surface $inc_rdzp -G${outfold}/rdzp_misfit_total_${region}.grd \
            -R$rdmin/$rdmax/$zpmin/$zpmax
    awk '{if ($3=='$rdbest'){print $1, $2, $4}}' ${outfold}/RMS_misfit_total_${region}.txt \
        | gmt surface $inc_Tpzp -G${outfold}/Tpzp_misfit_total_${region}.grd \
            -R$Tpmin/$Tpmax/$zpmin/$zpmax


    ref=1.05
    cont=$(gmt grdinfo ${outfold}/Tprd_misfit_total_${region}.grd -C | awk '{print $6*'$ref'}')
    gmt grdcontour ${outfold}/Tprd_misfit_total_${region}.grd -JX10c -R$Tpmin/$Tpmax/$rdmin/$rdmax -C+$cont \
        -D0.05_contour_Tprd_${region}.txt
    cont=$(gmt grdinfo ${outfold}/rdzp_misfit_total_${region}.grd -C | awk '{print $6*'$ref'}')
    gmt grdcontour ${outfold}/rdzp_misfit_total_${region}.grd -JX10c -R$rdmin/$rdmax/$zpmin/$zpmax -C+$cont \
        -D0.05_contour_rdzp_${region}.txt
    cont=$(gmt grdinfo ${outfold}/Tpzp_misfit_total_${region}.grd -C | awk '{print $6*'$ref'}')
    gmt grdcontour ${outfold}/Tpzp_misfit_total_${region}.grd -JX10c -R$Tpmin/$Tpmax/$zpmin/$zpmax -C+$cont \
        -D0.05_contour_Tpzp_${region}.txt

    awk 'NR>1{print $1, '$zpbest', $2, $3,"Tprd"}' 0.05_contour_Tprd_${region}.txt > ${outfold}/all_contours_${region}.txt
    awk 'NR>1{print '$Tpbest', $2, $1, $3,"rdzp"}' 0.05_contour_rdzp_${region}.txt >> ${outfold}/all_contours_${region}.txt
    awk 'NR>1{print $1, $2, '$rdbest', $3,"Tpzp"}' 0.05_contour_Tpzp_${region}.txt >> ${outfold}/all_contours_${region}.txt
    cont=$(awk 'NR==1 || $4<min {min=$4}END{print min}' ${outfold}/all_contours_${region}.txt)

    # Plot RMS ridge depth vs temperature:
    gmt grdcontour ${outfold}/Tprd_misfit_total_${region}.grd -R -J -C+$cont -D0.05_contour_Tprd_${region}.txt
    gmt makecpt -Cgray -T0.7/1.5/0.01 -I -D > misfit_depth.cpt
    gmt psbasemap -JX21/11.5c -R0/1/0/1 -B+t"${label_arr[$m]}" -K -Y20c > $ps
    gmt grdimage ${outfold}/Tprd_misfit_total_${region}.grd -JX10c -R$Tpmin/$Tpmax/$rdmin/$rdmax -Cmisfit_depth.cpt \
        -Bx100f${Tpinc}+l"Temperature (@~\260@~C)" -By200f${rdinc}+l"Ridge Depth (m)" -BWsNe -K -O >> $ps
    gmt grdcontour ${outfold}/Tprd_misfit_total_${region}.grd -R -J -C0.1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/Tprd_misfit_total_${region}.grd -R -J -A0.2 -Gn1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/Tprd_misfit_total_${region}.grd -R -J -C+$cont -K -O -Wthin,black,'-' >> $ps

    echo "1107.5 2030 @~\143@~@+2@+@-t@-" | gmt pstext -R -J -F+f12p+jLB -TO -Gwhite -W1p -O -K >> $ps
    awk '{if ($2=='$zpbest'){print $1, $3, $4}}' ${outfold}/RMS_misfit_total_${region}.txt | \
        awk 'NR==1 || $3<min {min=$3;Tp=$1;rd=$2}END{print Tp, rd}' | gmt psxy -R -J -Sc0.25c -Gred -Wthin -O -K >> $ps

    # Plot RMS ridge depth vs plate thickness:
    gmt grdcontour ${outfold}/rdzp_misfit_total_${region}.grd -R -J -C+$cont -D0.05_contour_rdzp_${region}.txt
    gmt grdimage ${outfold}/rdzp_misfit_total_${region}.grd -JX10c -R$rdmin/$rdmax/$zpmin/$zpmax -Cmisfit_depth.cpt -K -O \
        -Bx200f${rdinc}+l"Ridge Depth (m)" -By10f${zpinc}+l"Plate Thickness (km)" -BwsNE -X11c >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_total_${region}.grd -R -J -C0.1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_total_${region}.grd -R -J -A0.2 -Gn1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/rdzp_misfit_total_${region}.grd -R -J -C+$cont -K -O -Wthin,black,'-' >> $ps
    echo "2020 54 @~\143@~@+2@+@-t@-" | gmt pstext -R -J -F+f12p+jLB -TO -Gwhite -W1p -O -K >> $ps
    awk '{if ($1=='$Tpbest'){print $3, $2, $4}}' ${outfold}/RMS_misfit_total_${region}.txt | \
        awk 'NR==1 || $3<min {min=$3;rd=$1;zp=$2}END{print rd, zp}' | gmt psxy -R -J -Sc0.25c -Gred -Wthin -O -K >> $ps

    # Plot RMS plate thickness vs temperature:
    gmt grdinfo ${outfold}/Tpzp_misfit_total_${region}.grd -M
    gmt grdcontour ${outfold}/Tpzp_misfit_total_${region}.grd -R -J -C+$cont -D0.05_contour_Tpzp_${region}.txt
    gmt grdimage ${outfold}/Tpzp_misfit_total_${region}.grd -JX10c -R$Tpmin/$Tpmax/$zpmin/$zpmax -Cmisfit_depth.cpt -K -O \
        -Bx100f${rdinc}+l"Temperature (@~\260@~C)" -By10f${zpinc}+l"Plate Thickness (km)" -BWSnE -X-5.5c -Y-11c >> $ps
    gmt grdcontour ${outfold}/Tpzp_misfit_total_${region}.grd -R -J -C0.1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/Tpzp_misfit_total_${region}.grd -R -J -A0.2 -Gn1 -K -O -Wthin >> $ps
    gmt grdcontour ${outfold}/Tpzp_misfit_total_${region}.grd -R -J -C+$cont -K -O -Wthin,black,'-' >> $ps
    echo "1107.5 54 @~\143@~@+2@+@-t@-" | gmt pstext -R -J -F+f12p+jLB -TO -Gwhite -W1p -O -K >> $ps
    awk '{if ($3=='$rdbest'){print $1, $2, $4}}' ${outfold}/RMS_misfit_total_${region}.txt | \
        awk 'NR==1 || $3<min {min=$3;Tp=$1;zp=$2}END{print Tp, zp}' | gmt psxy -R -J -Sc0.25c -Gred -Wthin -O >> $ps

    gmt psconvert -P -Z $ps -Tj -E300 -A0.5c/0.5c



}


make_misfit_files(){

    echo -n > ${outfold}/grid_minima_${region}.txt
    gmt grdinfo ${outfold}/Tpzp_misfit_depth_${region}.grd -M -C | awk '{printf "%s %.0f %.0f\n", "depth", $12, $13}' > junk
    gmt grdinfo ${outfold}/rdzp_misfit_depth_${region}.grd -M -C | awk '{printf "%.0f %.5f\n", $12, $6}' > junk2
    paste junk junk2 >> ${outfold}/grid_minima_${region}.txt
    gmt grdinfo ${outfold}/Tpzp_misfit_heatflow_${region}.grd -M -C | awk '{printf "%s %.0f %.0f\n", "heatflow", $12, $13}' > junk
    gmt grdinfo ${outfold}/rdzp_misfit_heatflow_${region}.grd -M -C | awk '{printf "%.0f %.5f\n", $12, $6}' > junk2
    paste junk junk2 >> ${outfold}/grid_minima_${region}.txt
    gmt grdinfo ${outfold}/Tpzp_misfit_total_${region}.grd -M -C | awk '{printf "%s %.0f %.0f\n", "total", $12, $13}' > junk
    gmt grdinfo ${outfold}/rdzp_misfit_total_${region}.grd -M -C | awk '{printf "%.0f %.5f\n", $12, $6}' > junk2
    paste junk junk2 >> ${outfold}/grid_minima_${region}.txt

    gmt grdinfo ${outfold}/rdzp_misfit_depth_${Tp}_${region}.grd -M -C | \
        awk '{printf "%s %.0f %.0f %.0f %.5f\n", "depth_'${Tp}'", '$Tp', $13, $12, $6}' \
            >> ${outfold}/grid_minima_${region}.txt
    gmt grdinfo ${outfold}/rdzp_misfit_heatflow_${Tp}_${region}.grd -M -C | \
        awk '{printf "%s %.0f %.0f %.0f %.5f\n", "heatflow_'${Tp}'", '$Tp', $13, $12, $6}' \
            >> ${outfold}/grid_minima_${region}.txt
    gmt grdinfo ${outfold}/rdzp_misfit_total_${Tp}_${region}.grd -M -C | \
        awk '{printf "%s %.0f %.0f %.0f %.5f\n", "total_'${Tp}'", '$Tp', $13, $12, $6}' \
            >> ${outfold}/grid_minima_${region}.txt

}

make_grid(){

    inc_rough="-I0.2/1"
    inc_smooth="-I0.25/0.5"

    plate="${outfold}/Tgrid-${Tp}-${zp}-${rd}.dat"
    alpha=3.e-5
    cp=1.187e3
    g=9.81

    awk '{print $1, $3, $4}' $plate | gmt triangulate -Gtemp.grd $inc_smooth $rgn
    gmt grd2xyz temp.grd | awk '{if ($3=="NaN" && $2<50){print $1, $2, 0}else if \
        ($3=="NaN" && $2>=50) {print $1, $2, (('$Tp'+273)*exp(($2*1000*'$alpha'*'$g')/'$cp'))-273}else{print $0}}' \
            | gmt xyz2grd $inc_smooth -G$temperature_gmt $rgn
    rm -f temp.grd

}

temperature_shading () {

    # background temperature

    echo -e "0 0\n0 200\n200 200\n200 0" | gmt psxy $rgn $scale -G255/112.93/112.93 -K -O >> $ps

    # plot the temperature structure as shading

    gmt grdimage $temperature_gmt $rgn $scale -Ctemp.cpt -K -O >> $ps

}

temperature_labels () {

    gmt grdcontour $temperature_gmt $rgn $scale -C100 -L100/1500 -W0.5p,0 -A+f8p+n0c/0.16c -Gl0/120/100/0 -K -O >> $ps
    gmt grdcontour $temperature_gmt $rgn $scale -C+600 -Wthick -K -O >> $ps
    gmt grdcontour $temperature_gmt $rgn $scale -C+1000 -Wthick -K -O >> $ps

}

plot_all(){

	ps="${fold_plot_output}/temperature_output_${Tp}_${zp}_${rd}.ps"
	jpg="${fold_plot_output}/temperature_output_${Tp}_${zp}_${rd}.jpg"

	# Setup plot regins
	maxage=157
	rgnsub="-R0/${maxage}/0/8.5"
	rgnheat="-R0/${maxage}/2e-2/1e-0"
	scale_log="-JX12.25c/7.5cl"
	scalesub="-JX12.25c/-7.5c"
  rgnx="-R0/1/0/1"
  rgnx="-R-20/177/-1/9.5"
	scalex="-JX14.25c/9.5c"
    cptb="basement.cpt"
    gmt makecpt -T2/8/0.5 -Cpolar -D -I > $cptb

	# Plot subsidence trends:
	echo "subsidence"
    file_mod_subs="${outfold}/depth-${Tp}-${zp}-${rd}.dat"

	#Observations
    awk '{if ($6==0){print $4, $3, $3, $5}}' ${file_ref_subs} | \
        gmt psxy $scalesub $rgnsub -C$cptb -Ey2p/0.5p -Sc0.1c -Wthin -K -Y45c > $ps
    awk '{if ($6==1){print $4, $3, $3, $5}}' ${file_ref_subs} | \
        gmt psxy $scalesub $rgnsub -C$cptb -Ey2p/0.5p -St0.1c -Wthin -K -O >> $ps
    awk '{if ($6==2){print $4, $3, $3, $5}}' ${file_ref_subs} | \
        gmt psxy $scalesub $rgnsub -C$cptb -Ey2p/0.5p -Si0.1c -Wthin -K -O >> $ps

	#Models
	awk '{print $1, $2/1000}' $file_mod_subs | gmt psxy $rgnsub $scalesub -W1.5p,red -O -K >> $ps
	gmt psbasemap $rgnsub $scalesub -Bpxa20f10+l"Age (Ma)" -Bpya1f0.5+l"Water-Loaded Basement Depth (km)" -BWsNe -O -K >> $ps
  echo "3 0.3 a" | gmt pstext $rgnsub $scalesub -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps

	# Plot heat flow trends
	echo "heat flow"
	file_mod_hf="${outfold}/hf-${Tp}-${zp}-${rd}.dat"

	#Observations
    awk 'NR>1{printf "%.2f %.7f\n%.2f %.7f\n%.2f %.7f\n%.2f %.7f\n%s\n", \
        $1-1.25, $5, $1+1.25, $5, $1+1.25, $6, $1-1.25, $6, ">"}' $file_ref_hf | \
            gmt psxy $rgnheat $scale_log -Glightgray -W0.5p -X13.25c -O -K >> $ps
    awk 'NR>1{printf "%.2f %.7f\n%.2f %.7f\n%s\n", $1-1.25, $4, $1+1.25, $4, ">"}' $file_ref_hf | \
        gmt psxy $rgnheat $scale_log -W1p -O -K >> $ps

	#Models
	gmt psxy $file_mod_hf $rgnheat $scale_log -W1.5p,red -O -K >> $ps

	gmt psbasemap $rgnheat $scale_log -Bpx20f10+l"Age (Ma)" -By2f3+l"Heat Flow (W m@+-2@+)" -BwsNE -O -K >> $ps
  echo "3 0.3 b" | gmt pstext $rgnsub $scalesub -F+f18p+jTL -C+tO -Gwhite -W1p,black -N -O -K >> $ps

	#Plot isotherms HSC
	scale="-JX25.5c/-12.75c"
	rgn="-R0/${maxage}/0/158"

	temperature="${outfold}/Tgrid-${Tp}-${zp}-${rd}.grd"
	temperature_gmt=${outfold}/Tgrid-${Tp}-${zp}-${rd}.grd
    cpt="temp.cpt"

    gmt makecpt -T0/1400/100 -G-0.75/0.75 -Cpolar -D > $cpt

	gmt psbasemap $rgn $scale -Y-13.75c -X-13.25c -Ba0 -K -O >> $ps

    #rm -f $temperature_gmt
	if [ -e $temperature_gmt ]; then
	    echo "grid exists"
    else
        echo "grid doesn't exist, making..."
        make_grid
    fi

	temperature_shading
	temperature_labels

	gmt psbasemap $rgn $scale -Bxa20f10+l"Age (Ma)" -Bya20f5+l"Depth below Seafloor (km)" -BWeSn -O -K >> $ps
  echo "1.5 154.5 c" | gmt pstext $rgn $scale -F+f18p+jBL -C+tO -Gwhite -W1p,black -N -O >> $ps

	gmt psconvert $ps -Tj -A0.5c -P -Z


}

run_parameter_sweep(){
    # Run parameter sweep
    for Tp in $(seq $Tpmin $Tpinc $Tpmax); do

        for zp in $(seq $zpmin $zpinc $zpmax); do

            for rd in $(seq $rdmin $rdinc $rdmax); do

                # Run plate model
                plate_model

            done

        done

        # Check misfit
        assess_joint_misfit

    done
    #
    #
    ## Plot misfit results
    plot_all_misfit
    plot_joint_misfit
    #
}

run_mean_temperature_model(){
    ## Run model for mean mantle temperature
    Tp=$Tpcr

    for zp in $(seq $zpmin $zpinc $zpmax); do

        for rd in $(seq $rdmin $rdinc $rdmax); do

            # Run plate model
            continue
            plate_model

        done

    done

    # Check misfit for mean mantle temperature
    assess_joint_misfit_singleT

    # Plot misfit results for mean mantle temperature
    plot_singleT_misfit
}

make_output_best_fitting_model(){
    # Make temperature, depth and heatflow outputs for best-fitting models
    # Best-fit temperature model across all temperatures
    temp_flag=1
    Tp=$(awk '{if ($1=="total"){print $2}}' ${outfold}/grid_minima_${region}.txt)
    zp=$(awk '{if ($1=="total"){print $3}}' ${outfold}/grid_minima_${region}.txt)
    rd=$(awk '{if ($1=="total"){print $4}}' ${outfold}/grid_minima_${region}.txt)
    echo $Tp $zp $rd

    if [ -e ${codefold}/initial_geotherms/mac${Tp}.dat ]; then

        echo "Input temperature file already exists"

    else

        get_melting_profile

    fi

    # Make output files for best fit model across all temperatures
    plate_model

    # Plot output files
    plot_all
}

make_output_mean_temperature_model(){
    # Best-fit mean temperature models
    temp_flat=1
    Tp=1333
    zp=$(awk '{if ($1=="total_'${Tp}'"){print $3}}' ${outfold}/grid_minima_${region}.txt)
    rd=$(awk '{if ($1=="total_'${Tp}'"){print $4}}' ${outfold}/grid_minima_${region}.txt)

    if [ -e ${codefold}/initial_geotherms/mac${Tp}.dat ]; then

        echo "Input temperature file already exists"

    else

        get_melting_profile

    fi

    # Make output files for best fit model for mean temperature
    plate_model

    # Plot output files
    plot_all
}

copy_config(){

    cp config/jameshome.ini scripts/config.ini

}

####################################################################
# Set key parameters and directories
####################################################################

# Go to script fold
copy_config
cd scripts

# Set directories
fold_base=$(awk '$1 ~ /^base/' config.ini | awk '{print $3}')
fold_data_input=$(awk '$1 ~ /^data_input/' config.ini | awk '{print $3}')
fold_data_output=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')
fold_scripts=$(awk '$1 ~ /^scripts/' config.ini | awk '{print $3}')
fold_data_input=${fold_base}/${fold_data_input}
fold_data_output=${fold_base}/${fold_data_output}
fold_plot_output=${fold_base}/${fold_plot_output}
fold_scripts=${fold_base}/${fold_scripts}
input_heatflow=$(awk '$1 ~ /^input_heatflow/' config.ini | awk '{print $3}')
input_basement=$(awk '$1 ~ /^input_basement/' config.ini | awk '{print $3}')

codefold=${fold_scripts}/plate_model_code
meltfold=${fold_scripts}/melting_code
outfold=${fold_data_output}
hffold=${fold_data_input}/${input_heatflow}
subsfold=${fold_data_input}/${input_basement}
region="antarctic"
gmt gmtset FORMAT_FLOAT_OUT %.9f PS_MEDIA a0

# Set files
file_ref_subs=${subsfold}/${region}_basementdepth_1d.dat
file_ref_hf=${hffold}/2.5_Ma_binned_hf_stats_kappa_0.25_${region}_lt25Ma_edited.dat

####################################
# Set bounds for parameters sweep
####################################

# Zero-age ridge depth
rdmin=2000
rdmax=3000
rdinc=100

# Potential temperature
Tpmin=1100
Tpmax=1650
Tpinc=25

# Plate thickness
zpmin=70
zpmax=170
zpinc=10

####################################
# Set misfit parameters
####################################

# Subsidence data uncertainty
sigma_subs=700

# Data misfit weightings
weight_heatflow=1
weight_subs=1

############################
# Set additional parameters
############################

# Melting parameters
type=pd
rhoc=2890
rhocsd=40

# Mean mantle temperature
Tpcr=1333

############################
# Run the code
############################

# Make output folder
mkdir -p ${fold_data_output}

# Set temp flag to zero so large temperature file outputs are removed during parameter sweep
temp_flag=0

# run models
run_parameter_sweep
run_mean_temperature_model

# output misfit as a function of parameter sweep variables to find best-fitting and mean temperature optimal models
Tp=$Tpcr
make_misfit_files

# get best-fitting model data and plots
make_output_best_fitting_model
make_output_mean_temperature_model

# Clean up
rm -f input junk junk2 *.cpt resampled xpoints gmt.* 0.05_contour_*
