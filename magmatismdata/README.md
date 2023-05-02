This directory constructs a comparison between the estimated LAB depth structure and the geological record of past volcanism. 

data_input:
- conductive_1200C_isotherms are models of post-magmatic lithospheric rethickening (see Text S5 of Supporting Information)
- LAB_depth_models are various LAB depth models for comparison (Richards et al., 2020, JGR and Priestley et al., 2018, AGU)
- magmatism are records of past volcanism within the footprint of Antarctica from Sarbas et al. (2008, Geoinformatics 2008--Data to Knowledge) and Ball et al. (2021, Nature Communications)

data_output:
- filtered rock data bases based on age/location
- correlations between LAB depths and volcanism age records

plot_output:
- LAB depths, volcanism age records, and correlations between

scripts:
- extract_data.py filters the volcanism records
- make_LAB_magmatism_comparison.sh and calculate_spearman.py examines relationship between each LAB depth model and recorded ages of past volcanism, using MC sampling to incorporate uncertainty in each data set
- make_plot_LAB_magmatism.sh plots the results

References:
- Sarbas et al. (2008, Geoinformatics 2008--Data to Knowledge), link: https://gfzpublic.gfz-potsdam.de/pubman/item/item_10062 
- Ball et al. (2021, Nature Communications), DOI: 10.1038/s41467-021-22323-9
- Richards et al. (2020, JGR), DOI: https://doi.org/10.1029/2019JB019062
- Priestley et al. (2018, AGU), DOI: https://doi.org/10.1002/9781119249740.ch6