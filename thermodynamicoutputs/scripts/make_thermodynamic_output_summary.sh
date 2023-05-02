#!/bin/bash

get_comparison_data_output(){

    compare_Tp=$(awk '$1 ~ /^compare_potential_temperature/' config.ini | awk '{print $3}')
    compare_sol50=$(awk '$1 ~ /^compare_solidus_50km/' config.ini | awk '{print $3}')
    compare_date=$(awk '$1 ~ /^compare_date/' config.ini | awk '{print $3}')

    fold_actual_data_output=${fold_data_output_base}/comparison
    mkdir -p ${fold_actual_data_output}

    # make a new output dir if this combination of dates has not been looked at before
    # if it has, still make the plot but output it to the pre-exisiting dir (overwrite)
    output_digit=$(find ${fold_actual_data_output} -mindepth 1 -maxdepth 1 -type d | wc -l | awk '{print $1+1}' | awk '{printf "%.2d", $1}')
    outputs=$(find ${fold_actual_data_output} -mindepth 1 -maxdepth 1 -type d | sort)

    if [[ $(echo ${outputs} | awk '{print NF}' | awk '{if($1<1) print 999}') -eq 999 ]]; then
        
        fold_actual_data_output=${fold_actual_data_output}/comparison_${output_digit}
        mkdir -p ${fold_actual_data_output}
        f_log=${fold_actual_data_output}/comparison.log
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

                fold_actual_data_output=${output}

            else

                fold_actual_data_output=${fold_actual_data_output}/comparison_${output_digit}
                mkdir -p ${fold_actual_data_output}
                f_log=${fold_actual_data_output}/comparison.log
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

    echo ${fold_actual_data_output}

}

create_comparison_grids(){

    compare_Tp=$(awk '$1 ~ /^compare_potential_temperature/' config.ini | awk '{print $3}')
    compare_sol50=$(awk '$1 ~ /^compare_solidus_50km/' config.ini | awk '{print $3}')
    compare_date=$(awk '$1 ~ /^compare_date/' config.ini | awk '{print $3}')

    #fold_actual_data_output=${fold_data_output_base}/comparison
    mkdir -p ${fold_actual_data_output}

    fold_grid_in=${fold_data_output}/${grid_type}_grids/distribution/individual
    fold_data_output_comparison=${fold_data_output_base}/${compare_date}/Tp_${compare_Tp}_sol50_${compare_sol50}
    fold_grid_in_comparison=${fold_data_output_comparison}/${grid_type}_grids/distribution/individual

    fold_grid_out=${fold_actual_data_output}/${grid_type}_comparison_grids
    mkdir -p ${fold_grid_out}
    fold_grid_out=${fold_grid_out}/distribution
    mkdir -p ${fold_grid_out}
    fold_grid_out=${fold_grid_out}/individual
    mkdir -p ${fold_grid_out}

    idx_max=38

    for ((idx=38; idx<=$idx_max; idx++)); do

        path=${fold_grid_in}/depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_${idx}.grd
        path_comparison=${fold_grid_in_comparison}/depth_${depth}_Tp_${compare_Tp}_sol50_${compare_sol50}_idx_${idx}.grd
        path_out=${fold_grid_out}/depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_${idx}.grd
        if [[ ! -e $path_out ]]; then
            gmt grdmath $path $path_comparison SUB = $path_out
        else
            echo "exists!"
            gmt grdmath $path $path_comparison SUB = $path_out
        fi

        echo finished model ${idx} of ${idx_max}

    done

}

calculate_model_summary(){

    fold_grid_in=${fold_actual_data_output}/${grid_type}_grids/distribution/individual
    fold_grid_summary=${fold_actual_data_output}/${grid_type}_grids/distribution/summary
    mkdir -p ${fold_grid_summary}

    #if [[ -e ${fold_grid_summary}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd ]]; then
        #echo "file exists, skipping.."
        #return
    #fi

    rgn="-R0/360/-90/-60"
    inc="-I0.05d/0.05d"

    n_files=$(ls ${fold_grid_in}/*.grd | wc -l)
    idx_max=999

    idx=0
    echo mean $idx ${depth}
    path=${fold_grid_in}/depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_${idx}.grd

    gmt grdmath ${path} 0 DENAN = tmp_${grid_type}.grd
    gmt grdmath tmp_${grid_type}.grd 0 NEQ = tmp_data_count_${grid_type}.grd

    cp tmp_${grid_type}.grd sum_${grid_type}.grd
    cp tmp_data_count_${grid_type}.grd sum_data_count_${grid_type}.grd

    for ((idx=1; idx<=$idx_max; idx++)); do

        echo mean $idx ${depth}
        path=${fold_grid_in}/depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_${idx}.grd

        gmt grdmath ${path} 0 DENAN = tmp_${grid_type}.grd
        gmt grdmath tmp_${grid_type}.grd 0 NEQ = tmp_data_count_${grid_type}.grd

        gmt grdmath tmp_${grid_type}.grd sum_${grid_type}.grd ADD = sum_${grid_type}.grd
        gmt grdmath tmp_data_count_${grid_type}.grd sum_data_count_${grid_type}.grd ADD = sum_data_count_${grid_type}.grd

    done

    gmt grdmath sum_${grid_type}.grd sum_data_count_${grid_type}.grd DIV = mean_${grid_type}.grd

    idx=0
    echo std $idx ${depth}
    path=${fold_grid_in}/depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_${idx}.grd

    gmt grdmath ${path} 0 DENAN = tmp_${grid_type}.grd
    gmt grdmath tmp_${grid_type}.grd mean_${grid_type}.grd SUB SQR = dist_${grid_type}.grd

    for ((idx=1; idx<=$idx_max; idx++)); do

        echo std $idx ${depth}
        path=${fold_grid_in}/depth_${depth}_Tp_${Tp}_sol50_${sol50}_idx_${idx}.grd

        gmt grdmath ${path} 0 DENAN = tmp_${grid_type}.grd
        gmt grdmath tmp_${grid_type}.grd mean_${grid_type}.grd SUB SQR = tmp_dist_${grid_type}.grd
        gmt grdmath tmp_dist_${grid_type}.grd dist_${grid_type}.grd ADD = dist_${grid_type}.grd

    done

    gmt grdmath dist_${grid_type}.grd sum_data_count_${grid_type}.grd DIV SQRT = std_${grid_type}.grd

    cp sum_data_count_${grid_type}.grd ${fold_grid_summary}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_data_count.grd
    cp mean_${grid_type}.grd ${fold_grid_summary}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd
    cp std_${grid_type}.grd ${fold_grid_summary}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd
    #gmt grdfilter ${fold_grid_summary}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_mean.grd -D4 -Fg400 -G${fold_grid_summary}/${grid_type}_mean_filtered.grd
    #gmt grd2xyz ${fold_grid_summary}/${grid_type}_mean_filtered.grd > ${fold_grid_summary}/${grid_type}_mean_filtered.xyz
    #gmt grdfilter ${fold_grid_summary}/${grid_type}_depth_${depth}_Tp_${Tp}_sol50_${sol50}_std.grd -D4 -Fg400 -G${fold_grid_summary}/${grid_type}_std_filtered.grd

    rm tmp*.grd sum*.grd mean*.grd std*.grd dist*.grd

}

fold_data_output_base=$(awk '$1 ~ /^data_output/' config.ini | awk '{print $3}')
fold_plot_output=$(awk '$1 ~ /^plot_output/' config.ini | awk '{print $3}')
mkdir -p ${fold_data_output_base}
mkdir -p ${fold_plot_output}
Tp=$(awk '$1 ~ /^potential_temperature/' config.ini | awk '{print $3}')
sol50=$(awk '$1 ~ /^solidus_50km/' config.ini | awk '{print $3}')
date=$(awk '$1 ~ /^date/' config.ini | awk '{print $3}')

compare_Tp=$(awk '$1 ~ /^compare_potential_temperature/' config.ini | awk '{print $3}')
compare_sol50=$(awk '$1 ~ /^compare_solidus_50km/' config.ini | awk '{print $3}')
compare_date=$(awk '$1 ~ /^compare_date/' config.ini | awk '{print $3}')

mkdir -p ${fold_data_output_base}/${date}
mkdir -p ${fold_plot_output}/${date}
fold_data_output=${fold_data_output_base}/${date}/Tp_${Tp}_sol50_${sol50}
mkdir -p ${fold_data_output}
fold_actual_data_output=${fold_data_output} # unless changed for comparison calcs

for grid_type in temperature viscosity; do
    for depth in 75 150 250 350; do
        echo working on $grid_type grids at depth $depth km..
        #calculate_model_summary
        echo done
    done
done

for grid_type in temperature viscosity; do
    for depth in 75 150 250 350; do
        echo working on $grid_type grids at depth $depth km..
        #get_comparison_data_output
        #create_comparison_grids
        echo done
    done
done

for grid_type in temperature viscosity; do

    grid_type=${grid_type}_comparison

    for depth in 75 150 250 350; do
        echo working on $grid_type grids at depth $depth km..
        #get_comparison_data_output
        #calculate_model_summary
        echo done
    done
done