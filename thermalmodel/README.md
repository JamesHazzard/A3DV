This directory involves the preparation of a thermal model of Antarctic oceanic lithosphere (aka plate model), based on the formulation of Richards et al. (2018, JGR).

data_input:
- basement depth data
- heat flow data

data_output:
- results of grid-search optimization for best-fitting thermal parameters

plot_output:
- misfit and thermal model plots

scripts:
- melting_code scripts calculate melting profiles at a given potential temperature
- plate_model_code scripts calculate thermal model based on a given combination of temperature, ridge depth, and plate thickness
- main.sh runs a grid-search optimization for best-fitting thermal parameters by calling melting_code and plate_model_code scripts as necessary

references:
See Richards, F.D. et al. (2018, JGR), DOI: https://doi.org/10.1029/2018JB015998, and references therein