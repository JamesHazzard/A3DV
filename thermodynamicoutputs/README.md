This directory converts tomographic velocities into estimates of viscosity/density/temperature and calculates the mean/sdev of the resulting posterior distributions. 

data_input:
- anelastic_params are the outputs of the "posteriorsampling" directory, i.e. the *MAP* model and sieved parameter samples
- seismic_model are the depth slices of the tomographic model, i.e. ANT-20 from Lloyd et al. (2020, JGR)

data_output:
- viscosity/density/temperature grids

plot_output:
- plots of thermodynamic outputs at various depth slices

scripts:
- vs_to_thermo_conversions (see "parameteranalysis" directory) a la Richards et al. (2020, JGR)
- lib_quick_thermo.py for quick thermodynamic conversions
- libseis_min.py for reading of .grd data to python
- make_thermodynamic_output scripts used to propagate seismic velocity into thermodynamic variables using each anelasticity model, summarise the results, and plot the results

References:
- Lloyd, A.J. et al. (2020, JGR), DOI: https://doi.org/10.1029/2019JB017823 (ANT-20 tomography)
- Richards, F.D. et al. (2020, JGR), DOI: https://doi.org/10.1029/2019JB019062