This directory uses the calculated three-dimensional temperature structure to calculate spatially variable geothermal structure beneath Antarctica, and from this estimate lithosphere-asthenosphere boundary (LAB) depth and geothermal heat flow (GHF).

data_input:
- anelastic_params are outputs of "posteriorsampling" directory, i.e. a summary of inverted parameters
- crustal_grids are basement depth / crustal thickness / Moho discontinuity depth grids, consistent with the tomography model (Lloyd et al., 2020, JGR). Oceanic age grid is also here (Richards et al., 2018, JGR).
- LAB1200_to_baseTBL is a calibrated conversion between LAB depth (according to the 1200 degree C isotherm) and the base of the thermal boundary layer (TBL)
- sample_locations is a set of locations we wish to calculate and plot geotherms at
- seismic_model is the downsampled tomography model (Lloyd et al., 2020, JGR)

data_output:
- geotherms are calculated temperature-depth structures
- LAB1200 contains spatially variable LAB depth, according to the 1200 degree C isotherm
- HF_const_Tp and baseTBL_const_Tp are geotherm fitting outputs calculated while ignoring any spatial variation in mantle potential temperature (see manuscript for details)
- LAB1200_to_baseTBL calculates the empirically fitted relationship between the 1200 degree C isotherm and the baseTBL_const_Tp 
- HF are geotherm fitting outputs calculated while accounting for spatial variation in mantle potential temperature, which requires the above empirical conversion between LAB depth and base TBL depth as an input (hence its appearance in the data_input directory)

scripts:
- fit_geotherm is a Fortran 90 script allowing us to fit a geotherm to the seismically inferred temperature structure (Mckenzie et al., 2005, EPSL)
- lib_anelasticity.py (with libseis_min.py) allows reading of a tomography model and conversion into an array of temperatures (Richards et al., 2020, JGR)
- make_output_dir.py adds an output folder according to the date/time of analysis run
- make_geotherms.py allows fitting of geotherms at specific locations, ready for plotting
- make_LAB1200.py allows measurement of spatially variable depth to the 1200 degree C isotherm
- make_LAB1200_to_baseTBL_correlation.py calculates empirical relationship between LAB1200 and baseTBL depths
- make_HF.py allows geotherm fitting and thus estimate of heat flow into the base of the Antarctic Ice Sheet
- submit_hpc_HF.sh is an example of a supercomputer job submission script that might be used to send geotherm fitting jobs off to a cluster

References:
- Lloyd et al. (2020, JGR), DOI: https://doi.org/10.1029/2019JB017823
- Richards et al. (2018, JGR), DOI: https://doi.org/10.1029/2018JB015998
- McKenzie et al. (2005, EPSL), DOI: https://doi.org/10.1016/j.epsl.2005.02.005
- Richards et al. (2020, JGR), DOI: https://doi.org/10.1029/2019JB019062