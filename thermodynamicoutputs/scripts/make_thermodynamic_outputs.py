from libseis_min import load_grd_file, lay_seis_model
from lib_quick_thermo import anelasticity_model
import numpy as np
from glob import glob
from os import sys
from os.path import join
import configparser
import os
import time
import sys, getopt
from datetime import datetime
import xarray as xr

model = 'ANT-20' # update if using a different seismic model
model_data = 'ant_region_25-400km_0.05d_60S'
t_start = time.time()

###############################################################################################
# 1. GET DIRECTORIES FROM CONFIG
###############################################################################################

config_obj = configparser.ConfigParser()
config_obj.read('config.ini')
variables = config_obj["variables"]
potential_temperature = variables["potential_temperature"]
solidus_50km = int(variables["solidus_50km"])
date = variables["date"]
folds = config_obj["directories"]
fold_base = folds["base"]
fold_data_input = folds["data_input"]
fold_anelastic_params = os.path.join(fold_data_input, folds["input_anelastic_params"], date, f'Tp_{potential_temperature}_sol50_{solidus_50km}')
fold_seismic_model = os.path.join(fold_data_input, folds["input_seismic_model"], model_data)

fold_data_output = os.path.join(folds["data_output"], date, f'Tp_{potential_temperature}_sol50_{solidus_50km}')
data_output_types = ['density', 'temperature', 'viscosity']

os.makedirs(fold_data_output, exist_ok=True)
for o_type in data_output_types:

    os.makedirs(os.path.join(fold_data_output, f'{o_type}_grids'), exist_ok=True)
    os.makedirs(os.path.join(fold_data_output, f'{o_type}_grids', 'distribution'), exist_ok=True)
    os.makedirs(os.path.join(fold_data_output, f'{o_type}_grids', 'distribution', 'individual'), exist_ok=True)
    os.makedirs(os.path.join(fold_data_output, f'{o_type}_grids', 'distribution', 'summary'), exist_ok=True)
    os.makedirs(os.path.join(fold_data_output, f'{o_type}_grids', 'MAP'), exist_ok=True)

###############################################################################################
# 2. IMPORT ANELASTICITY PARAMETERS BASED ON USER INPUT
###############################################################################################

# initiate variables
fitfile_type = ''
fitfile_index = -1

# read -v variable for type (specify "MAP" for MAP) and -i variable for index (if using distribution)
opts, args = getopt.getopt(sys.argv[1:], "v:i:")
for opt,arg in opts:
    if opt == '-v':
        fitfile_type = str(arg)
    if opt == '-i':
        fitfile_index = eval(arg)

if fitfile_type.strip().lower() in ['map', 'maximum a posteriori', 'most probable']:

    fitfile = os.path.join(fold_anelastic_params, 'MAP_model.txt')
    pars=np.loadtxt(fitfile)[1:]
    print('using anelasticity parameters from MAP model..')
    print('parameters are', pars)
    fitfile_type = 'MAP'
    fitfile_index = 'N/A'

elif fitfile_index >= 0:

    fitfile = os.path.join(fold_anelastic_params, 'data_sample.txt')
    pars=(np.loadtxt(fitfile).T)[1:]
    pars=pars[:,fitfile_index]
    print('using anelasticity parameters from sample', fitfile_index, 'of distribution..')
    fitfile_type = 'distribution'
    fitfile_index = str(fitfile_index)

###############################################################################################
# 3. IMPORT SEISMIC MODEL
###############################################################################################

# Initialise conversion
main_inversion  = anelasticity_model(eta_0=float(10**pars[3]), mu0=float(pars[0]*1.e9), dmudT=float(pars[1]*1.e9), dmudp=float(pars[2]),\
                  act_eng=float(pars[4]), act_vol=float(pars[5]), solgrad=float(pars[6]), sol50=solidus_50km, alpha_B=0.38,\
                                    A_B=0.664, delphi=0.0, temps=np.linspace(273.,4473.,4201), dens_flg=1);
# Read in seismic file
seis_pth = fold_seismic_model
# Get the names of all the files
seis_files = glob(join(seis_pth,'*.grd'))

# Make an array consisting of all the depths
seismic_depths = np.zeros(len(seis_files));
# Set maximum depth
maxdep=400.
# Obtain the depths of the seismic models
for i in range(len(seis_files)):
   seismic_depths[i] = float(seis_files[i][seis_files[i].rfind('/')+1:seis_files[i].rfind('km')]);

# Sort the layers and remove slices deeper than maxdep (400 km)
seismic_depths.sort(); seismic_depths = seismic_depths[seismic_depths<=maxdep]*1.e3;

# Define an array for seismic layers
seismic_layers_orig    = [None]*len(seismic_depths);
seismic_layers    = [None]*len(seismic_depths);
seis_object    = [None]*len(seismic_depths);

# Go through all the seismic_depths, and read in the files
for i in range(len(seismic_depths)):
   print(str("Reading in depth %4.0f km" %(seismic_depths[i]/1.e3)))
   seis_object[i] = lay_seis_model(model_path = join(seis_pth,\
               str("%ikm.grd" %(seismic_depths[i]/1.e3))));

   # Write data
   seismic_layers_orig[i] = seis_object[i].data;
   seismic_layers[i] = seis_object[i].data;

   # Check for NaNs and set to layer mean
   seismic_layers[i][seismic_layers_orig[i]!=seismic_layers_orig[i]]=np.nanmean(seismic_layers_orig[i]);

# Get lon/lat for depth slice
lat = seis_object[i].lat
lon = seis_object[i].lon
n_lat = len(lat)
n_lon = len(lon)
n_locs = int(len(lat) * len(lon))
locs = np.stack(np.meshgrid(lon, lat)).reshape(2, n_locs)

###############################################################################################
# 5. CONVERT SEISMIC VELOCITY INTO TEMPERATURE, DENSITY AND VISCOSITY
###############################################################################################

def T2Tp(z,T):
    alpha=3.e-5
    Cp=1187.
    g=9.81
    Tp=(T+273.15)*np.exp(-(alpha*g*z*1000.)/Cp)-273.15
    return np.array([z,Tp])

# Convert the seismic model to temperature
temp_fld, dens_fld, visc_fld = main_inversion.vs_2_thermo(depths=seismic_depths,\
                                    seis_model=seismic_layers);

###############################################################################################
# 6. FLATTEN ARRAYS AND SAVE
###############################################################################################

data_output_types = ['density', 'temperature', 'viscosity']
data_output_types_shorthand = ['dens', 'temp', 'visc']

for i in range(len(seismic_depths)):

    depth_str = str(int(seismic_depths[i]/1.e3))

    tmp_data = temp_fld[i] - 273
    da_temp = xr.DataArray(
        data = tmp_data,
        dims = ['Latitude', 'Longitude'],
        coords = {'Latitude': lat, 'Longitude': lon},
        name = data_output_types[1],
        attrs = {'actual_range': [np.min(tmp_data), np.max(tmp_data)]}
    )

    tmp_data = dens_fld[i]
    da_dens = xr.DataArray(
        data = tmp_data,
        dims = ['Latitude', 'Longitude'],
        coords = {'Latitude': lat, 'Longitude': lon},
        name = data_output_types[0],
        attrs = {'actual_range': [np.min(tmp_data), np.max(tmp_data)]}
    )

    tmp_data = np.log10(visc_fld[i])
    da_visc = xr.DataArray(
        data = tmp_data,
        dims = ['Latitude', 'Longitude'],
        coords = {'Latitude': lat, 'Longitude': lon},
        name = data_output_types[2],
        attrs = {'actual_range': [np.min(tmp_data), np.max(tmp_data)]}
    )

    output_xarrays = [da_dens, da_temp, da_visc]

    #for idx_output in range(len(data_output_types)):
    for idx_output in [1, 2]:

        o_type = data_output_types[idx_output]
        print(f"saving {o_type}..")
        o_file_name = os.path.join(fold_data_output, f'{o_type}_grids', 'distribution', 'individual',\
                                     f'depth_{depth_str}_Tp_{potential_temperature}_sol50_{solidus_50km}_idx_{fitfile_index}.grd')
        xarray = output_xarrays[idx_output]
        xarray.to_netcdf(o_file_name, encoding={o_type: {'dtype': 'float32', 'zlib': True}})                            

t_end = time.time()
print(t_end - t_start, "seconds from start to finish")