This directory allows us to analyse and plot the distribution of posterior anelasticity models. 

Data input is taken from BANCAL22 run output directory. Use the config file to adjust where you wish to pull inversion output data from. 

data_output:
- credible intervals showing fit between inversion data and posterior models
- histograms showing relationship between *a priori* and *a posteriori* parameters
- trade_offs describing the covariance between paris of parameters

plot_output:
- plots allowing us to visualise the data outputs

scripts:
- vs_to_thermo_conversions contains Fortran 90 scripts needed to convert seismic velocity into thermodynamic parameters using the YT16 (Yamauchi and Takei, 2016, JGR) anelasticity model. These scripts were developed in association with the work of Richards et al. (2020, JGR)
- make_credible_interval... scripts calculate and plot the relationship between inversion data and model
- make_parameter_comparison.py allows comparison of parameters between inversion runs
- make_parameter_histograms.sh allows plot of parameter histograms
- make_parameter_summary.py is used to calculate simple summaries of the parameter outputs
- make_trade_off_density_plot.sh allows calculation and plotting of sampling density between any given two parameters within the anelasticity parameterisation, which serves as a statistical proxy for covariance
- main.sh weaves scripts together and allows control of overall directory activity

References:
- Yamauchi and Takei (2016, JGR), DOI: https://doi.org/10.1002/2016JB013316
- Richards, F.D. et al. (2020, JGR), DOI: https://doi.org/10.1029/2019JB019062