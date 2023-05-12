# A3DV (Antarctica 3D Viscosity)

The code provided in this repository relates to the manuscript entitled "Probabilistic Assessment of Antarctic Thermomechanical Structure: Impacts on Ice Sheet Stability". Something not clear about how to use this code? Interested in working with me? Email me: **j.hazzard20@imperial.ac.uk**. 

Please respect the user license associated with this code. You should cite the associated manuscript (DOI: https://doi.org/10.1029/2023JB026653) in any research produced using the materials provided herein.

## Software languages

The code in this directory is mostly bash script and Python, with a small component of Fortran 90. All code was written and optimised for use with my Ubuntu OS, and Python 3.10. Small adjustments may be necessary to make the code work with your system. 

## The workflow 

This guide explains how to use the code provided in this repository in order to generate the outputs described in the manuscript. Let's start with a summary of what can be found here. 

### Code directories
1. **thermalmodel** (Preparation of thermal model of Antarctic oceanic lithosphere aka plate model)
2. **inversiondata** (Preparation of Bayesian inversion data constraints, namely: plate, adiabat, attenuation and viscosity)
3. **BANCAL22** (Bayesian inversion procedure to calibrate experimental parameterisations of anelasticity)
4. **posteriorsampling** (Sieving of posterior distribution of anelasticity parameters)
5. **parameteranalysis** (Analysis and plotting of distribution of anelasticity parameters: credible interval data fits, histograms, trade-offs)
6. **thermodynamicoutputs** (Conversion of Vs into estimates of viscosity/density/temperature and calculate mean/sdev of resulting distributions)
7. **geothermfitting** (Fitting of 3D temperature structure to estimate lithosphere-asthenosphere boundary (LAB) depth and geothermal heat flow (GHF))
8. **magmatismdata** (Comparison of estimated LAB structure to geological record of past volcanism)
9. **pappa19** (Comparison of estimated LAB structure to previous study)
10. **westvseast** (Separation of calculated thermomechanical outputs into West vs. East Antarctica for analysis)
11. **complexviscosity** (Analysis of the effect of time-dependence on agreement between seismic- and geodetic-based inferences of Antarctic viscosity)
### Data directories
1. **ANT-20_S_i_abs** (ANT-20 seismic data)
2. **QRFSI12** (QRFSI12 attenuation data)
3. **BANCAL22_runs** (inversion output data)

Code directories contain scripts that turn an input into an output. Data directories store data, which are used at various points within the code.

### Code directory structure

Each code directory contains a **config** folder, in which you can store a [config_name].ini file. An example of how the config file should be structured is included in each code directory. The config file tells the scripts in the code directory where it can find inputs, and where it should store outputs. 

You will also find **data_input**, **data_output**, **plot_output** and **scripts** folders. These are relatively self-explanatory. The names of these folders match the example config file.

Finally, each code directory contains a **main.sh** file, which is a Bash script you can use to run all of the functionality included in the directory of interest.