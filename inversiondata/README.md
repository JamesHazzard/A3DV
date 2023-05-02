This directory involves the preparation of Bayesian inversion data constraints: plate, adiabat, attenuation, and viscosity. 

data_input:
- exclusion polygons to mask out certain areas of oceanic crust (Hoggard et al., 2017, JGR)
- ocean age grids (Richards et al., 2018, JGR)
- plate thermal model grid (Richards et al., 2018, JGR - see "thermalmodel" directory)

data_output:
- prepared inversion data

scripts:
- make_adiabat.sh ties tomographic velocities to oceanic age constraints and constructs an adiabatic Vs-T data set
- make_plate.sh uses pre-determined (see above) velocity-age-depth constraint and co-variation with plate thermal model temperature to construct a plate Vs-T data set
- make_attenuation.sh ties seismic attenuation to tomographic velocity to construct an attenuation data set
- make_viscosity.sh constructs a bulk shear-viscosity estimate data set
- make_BANCAL22_data.sh prepares above data sets in BANCAL22 inversion-ready format and adds information relating to the solidus and potential temperatures to use in the inversion. It can also transfer the prepared data into the BANCAL22 inversion directory of choice, ready for use

references:
- Hoggard, M.J. et al. (2017, JGR), DOI: https://doi.org/10.1002/2016JB013457
- Richards, F.D. et al. (2018, JGR), DOI: https://doi.org/10.1029/2018JB015998