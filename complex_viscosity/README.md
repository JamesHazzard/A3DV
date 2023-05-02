This directory allows for analysis of the effect of time-dependence on the agreement between seismically and geodetically inferred Antarctic upper mantle viscosity (Barletta et al., 2018, Samrat et al., 2020, Ivins et al., 2011, Wolstencroft et al., 2015).

data_input:
- anelastic_params is the set of posterior anelasticity models from the "posteriorsampling" directory output
- sample_locs is a set of locations where we wish to sample viscosity
- temperature_grids contains seismically inferred temperatures at the study locations, as sampled from the outputs of "thermodynamicoutputs" directory (this input is created by the script make_sample_loc_grids.sh)

data_output:
- apparent_viscosity records the "apparent" or "effective" viscosity relevant to the combination of seismically inferred steady-state viscosities, the transient rheology YT16 (Yamauchi and Takei, 2016), and the ice history assumed in the GPS study we are comparing our results to
- strain_history records the associated synthetic strain rate observations 

plot_output:
- viscosity_kdes shows location of samples and agreement between steady-state, time-dependent and geodetically inferred viscosities

scripts:
- make_loading_history.py models the ice loading history of the GPS studies used in the analysis
- make_relaxation_spectrum.py models the relationship between steady-state and time-dependent viscosity, which depends on the modelled loading history (based on the formulation of Lau et al., 2021)
- make_sample_loc_grids.sh constructs the necessary input temperature grids
- make_sample_loc_plot.sh plots the locs of the GPS studies
- make_viscosity_summary.py summarises the outputs of the viscosity analysis and produces a kernel density estimate (kde)
- make_viscosity_summary_plots.sh produces plots showing the comparison between seismic and geodetic inferences of viscosity across spatial and temporal scales

References:
- Lau et al. (2021, JGR), DOI: 10.1029/2021JB022622
- Barletta et al. (2018, Science), DOI: 10.1126/science.aao1447
- Samrat et al. (2020, GJI), DOI: 10.1093/gji/ggaa229
- Ivins et al. (2021, JGR), DOI: https://doi.org/10.1029/2010JB007607
- Wolstencroft et al. (2015, GJI), DOI: 10.1093/gji/ggv327