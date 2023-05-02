from libseis_min import load_grd_file, lay_seis_model
from lib_anelasticity import anelasticity_model
import numpy as np
import scipy
from glob import glob
from os import sys
from os.path import join
import configparser
import os
import subprocess
import time
import sys, getopt
from datetime import datetime

model = 'ANT-20' # update if using a different seismic model
model_data = 'ant_cont_region_25-400km_downsampled_0.5d'

###############################################################################################
# 1. GET DIRECTORIES FROM CONFIG
###############################################################################################

config_obj = configparser.ConfigParser()
config_obj.read('config.ini')
variables = config_obj["variables"]
potential_temperature = variables["potential_temperature"]
solidus_50km = int(variables["solidus_50km"])
folds = config_obj["directories"]
fold_base = folds["base"]
fold_data_input = folds["data_input"]
fold_anelastic_params = os.path.join(fold_data_input, folds["input_anelastic_params"], f'Tp_{potential_temperature}_sol50_{solidus_50km}')
fold_crustal_grids = os.path.join(fold_data_input, folds["input_crustal_grids"])
fold_seismic_model = os.path.join(fold_data_input, folds["input_seismic_model"], model_data)
fold_LAB_to_baseTBL = os.path.join(fold_data_input, folds["output_LAB1200_TBL"])
fold_locs = os.path.join(fold_data_input, folds["input_locs"])

now_date = datetime.strptime(variables["date"], "%Y-%m-%d_%H-%M-%S")
now_date_string = f'RUN_{now_date.strftime("%Y-%m-%d_%H-%M-%S")}'
year = f'year_{now_date.strftime("%Y")}'
month = f'month_{now_date.strftime("%m")}'
day = f'day_{now_date.strftime("%d")}'

fold_data_output = os.path.join(folds["data_output"], f'Tp_{potential_temperature}_sol50_{solidus_50km}', 'archive', now_date_string)
fold_geotherms = os.path.join(fold_data_output, folds["output_geotherms"])
fold_geotherms_latest = os.path.join(folds["data_output"], f'Tp_{potential_temperature}_sol50_{solidus_50km}', 'latest', folds["output_geotherms"])

os.makedirs(fold_geotherms, exist_ok=True)
os.makedirs(fold_geotherms_latest, exist_ok=True)

###############################################################################################
# 2. IMPORT EMPIRICAL CONVERSION BETWEEN LAB1200 AND BASETBL TO ESTIMATE POTENTIAL TEMPERATURE
###############################################################################################

file_LAB_to_baseTBL = os.path.join(fold_LAB_to_baseTBL, f'LAB1200_to_baseTBL_conversion_Tp_{potential_temperature}_sol50_{solidus_50km}.txt')
conv_data = np.loadtxt(file_LAB_to_baseTBL)

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
    outfold = os.path.join(fold_geotherms, 'MAP')
    outfold_latest = os.path.join(fold_geotherms_latest, 'MAP')
    os.makedirs(outfold, exist_ok=True)
    os.makedirs(outfold_latest, exist_ok=True)
    pars=np.loadtxt(fitfile)[1:]
    print('using anelasticity parameters from MAP model..')
    print('parameters are', pars)
    fitfile_type = 'MAP'
    fitfile_index = 'N/A'

elif fitfile_index >= 0:

    fitfile = os.path.join(fold_anelastic_params, 'data_sample.txt')
    outfold = os.path.join(fold_geotherms, 'distribution')
    outfold_latest = os.path.join(fold_geotherms_latest, 'distribution')
    os.makedirs(outfold, exist_ok=True)
    os.makedirs(outfold_latest, exist_ok=True)
    outfold = os.path.join(outfold, 'individual')
    outfold_latest = os.path.join(outfold_latest, 'individual')
    os.makedirs(outfold, exist_ok=True)
    os.makedirs(outfold_latest, exist_ok=True)
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

###############################################################################################
# 4. IMPORT CRUSTAL GRIDS
###############################################################################################

# Load the age grid map, where continents are set to "NaN"
cnt_map_lat, cnt_map_lon, cnt_map_data = \
            load_grd_file(os.path.join(fold_crustal_grids, 'ant_cont_region_sampled_oc_age_antarctica_0.5d.grd'))

# Load crustal thickness map
crst_map_lat, crst_map_lon, crst_map_data = load_grd_file(os.path.join(fold_crustal_grids, 'ant_cont_region_sampled_modified_crust1_crsthk_0.5d.grd'))

# Load moho depth map
moho_map_lat, moho_map_lon, moho_map_data = load_grd_file(os.path.join(fold_crustal_grids, 'ant_cont_region_sampled_modified_crust1_moho_0.5d.grd'))

# Load basement depth map
basement_map_lat, basement_map_lon, basement_map_data = load_grd_file(os.path.join(fold_crustal_grids, 'ant_cont_region_sampled_modified_crust1_basement_0.5d.grd'))

###############################################################################################
# 5. CONVERT SEISMIC VELOCITY INTO TEMPERATURE
###############################################################################################

def T2Tp(z,T):
    alpha=3.e-5
    Cp=1187.
    g=9.81
    Tp=(T+273.15)*np.exp(-(alpha*g*z*1000.)/Cp)-273.15
    return np.array([z,Tp])

# Convert the seismic model to temperature
temp_fld = main_inversion.vs_2_temp(depths=seismic_depths,\
                                    seis_model=seismic_layers);

# Make temperature cube: dimension 1: depths,; 2: lon; 3: lat
temp_cube=np.zeros((len(seismic_depths),np.shape(temp_fld[0])[0],np.shape(temp_fld[0])[1]))
for i in range(len(seismic_depths)):
    temp_cube[i,:,:]=temp_fld[i]

###############################################################################################
# 6. SET UP ARRAYS NEEDED TO STORE INPUTS & OUTPUTS OF LOCATION LOOP
###############################################################################################

# Make location and initial temperature array
ndepths=len(seismic_depths)
locs=np.zeros((np.size(temp_fld[0]),2))
temp_out=np.zeros((np.size(temp_fld[0]),ndepths))
k=0
for i in range(len(cnt_map_lat)):
    for j in range(len(cnt_map_lon)):
        locs[k,0]=cnt_map_lon[j]
        locs[k,1]=cnt_map_lat[i]
        for zo in range(ndepths):
            temp_out[k,zo]=temp_fld[zo][i,j]
        k=k+1

# Make crustal thickness array
k=0
nlocs=np.size(temp_fld[0])
crust_arr=np.zeros(nlocs)
moho_arr=np.zeros(nlocs)
basement_arr=np.zeros(nlocs)
for i in range(len(crst_map_lat)):
    for j in range(len(crst_map_lon)):
        crust_arr[k]=crst_map_data[i,j]
        moho_arr[k]=moho_map_data[i,j]
        basement_arr[k]=basement_map_data[i,j]
        k=k+1

# Make depth and temperature array
nlocs=np.size(temp_fld[0])
temp_arr=np.zeros((2,ndepths+1,nlocs))
for t in range(np.shape(locs)[0]):
    temp_arr[0,1:,t]=seismic_depths[:]*1e-3
    temp_arr[1,1:,t]=temp_out[t,:]-273.15

# Set crustal thickness,heat production & conductivity parameters
# Set conductivity parameters
kmod=np.full(np.shape(cnt_map_data),8)
kmod[cnt_map_data==cnt_map_data]=8
kmod=np.reshape(kmod,np.size(cnt_map_data))
# Set crustal conductivity
crustk=np.full(np.shape(cnt_map_data),2.5)
crustk[cnt_map_data==cnt_map_data]=2.59
crustk=np.reshape(crustk,np.size(cnt_map_data))
# Set crustal heat production
h_upper=np.full(np.shape(cnt_map_data),1.0e-6)
h_upper[cnt_map_data==cnt_map_data]=0.
h_upper=np.reshape(h_upper,np.size(cnt_map_data))
h_lower=np.full(np.shape(cnt_map_data),0.3e-6)
h_lower[cnt_map_data==cnt_map_data]=0.
h_lower=np.reshape(h_lower,np.size(cnt_map_data))
# Set viscosity parameters
vis=np.full(np.size(cnt_map_data),9.e16)
# Set potential temperature parameters (overwritten by regional estimate)
tp=np.full(np.size(cnt_map_data),float(potential_temperature)) # in const Tp case should use same Tp as where params came from
# Set crustal thickness parameters
t=crst_map_data
t=np.reshape(crst_map_data,np.size(cnt_map_data))

# Set MBL array
MBL=np.zeros(np.shape(temp_arr)[2])
# Set TBL array
TBL=np.zeros(np.shape(temp_arr)[2])
# Set LAB array
LAB=np.zeros(np.shape(temp_arr)[2])
# Set LAB1200 array
LAB1200=np.zeros(np.shape(temp_arr)[2])
# Set estimated base TBL array
estimate_baseTBL=np.zeros(np.shape(temp_arr)[2])
# Set base TBL array
baseTBL=np.zeros(np.shape(temp_arr)[2])
# Set misfit array
misfit=np.zeros(np.shape(temp_arr)[2])
# Set Moho temperature
Tmoho=np.zeros(np.shape(temp_arr)[2])
# Set surface heat flux array
Sflux=np.zeros(np.shape(temp_arr)[2])
# Set Moho heat flux array
Mflux=np.zeros(np.shape(temp_arr)[2])
# Remove old geotherm output file
subprocess.Popen('rm geoth.out', shell=True,stdout=subprocess.PIPE)

# Set depth increment for interpolation
zinc=1.

# Set LAB temperature
TLAB=1200.

t_start = time.time()

###############################################################################################
# 7. EXTRACT LOCATION INDICES WE WISH TO SAMPLE
###############################################################################################

data_locs = np.loadtxt(os.path.join(fold_locs, "sample_locs.txt"))
n_locs = np.shape(data_locs)[0]
idx_sample_loc = []

for idx_loc in range(n_locs):

    lon = data_locs[idx_loc,0]
    lat = data_locs[idx_loc,1]
    idx_sample_loc.append(np.where((locs[:,0] == lon) & (locs[:,1] == lat))[0].item())

###############################################################################################
# 8. LOOP OVER LOCATIONS IN GRID AND CALCULATE LAB1200 AT EACH POINT
###############################################################################################

def main():

    k = 1

    for zl in idx_sample_loc:

        print(f"working on geotherm {k} of {len(idx_sample_loc)}")

        # Interpolate temperature profile
        lon = locs[zl,0]
        lat = locs[zl,1]
        if os.path.isfile(os.path.join(outfold_latest, f'lon_{lon}_lat_{lat}_geoth_6.txt')):
            continue
        np.savetxt(os.path.join(outfold_latest, f'lon_{lon}_lat_{lat}_geoth_1.txt'), temp_arr[:,:,zl].T, fmt='%3.0f %12.5f')
        mantle_arr = temp_arr[:,np.where(temp_arr[0,:,zl]>=moho_arr[zl])[0],zl] # Get rid of any data from the crust
        np.savetxt(os.path.join(outfold_latest, f'lon_{lon}_lat_{lat}_geoth_2.txt'), mantle_arr.T, fmt='%3.0f %12.5f')
        mantle_arr = mantle_arr[:,np.where(mantle_arr[0,:]>=mantle_arr[0,np.gradient(mantle_arr[:,:])[1][1]>=0][0])][:,0] # Get rid of any data shallower than an inflection point
        np.savetxt(os.path.join(outfold_latest, f'lon_{lon}_lat_{lat}_geoth_3.txt'), mantle_arr.T, fmt='%3.0f %12.5f')
        surface_arr_depths = np.array([basement_arr[zl]])
        surface_arr_T = np.array([0.0])
        surface_arr = np.stack((surface_arr_depths, surface_arr_T))
        arr = np.concatenate((surface_arr, mantle_arr), axis = 1)
        np.savetxt(os.path.join(outfold_latest, f'lon_{lon}_lat_{lat}_geoth_4.txt'), arr.T, fmt='%12.5f %12.5f')
        minz = arr[0,0]
        maxz = arr[0,-1]
        trial_depths = np.linspace(minz, maxz, int((maxz-minz)/zinc)+1)
        intp=scipy.interpolate.Akima1DInterpolator(arr[0,:],arr[1,:]) # [depth/T ; depth slice] Interpolate remaining profile
        arr=np.stack((trial_depths,intp(trial_depths)))
        np.savetxt(os.path.join(outfold_latest, f'lon_{lon}_lat_{lat}_geoth_5.txt'), arr.T, fmt='%12.5f %12.5f')

        arr_LAB = arr.copy().T
        if np.size(arr_LAB[np.where((arr_LAB[:,1]>=TLAB)),:]) > 0:
            LAB1200[zl] = arr_LAB[np.where((arr_LAB[:,1]>=TLAB)),:][0][0,0]
            estimate_baseTBL[zl] = (LAB1200[zl] * conv_data[0]) + conv_data[1]
            find = np.where((arr[0] >= estimate_baseTBL[zl]))[0]
            if np.size(find) < 1:
                belowTBL = arr[:,-1]
                belowTBL_Tp = T2Tp(belowTBL[0], belowTBL[1])
            elif np.size(find) >= 1:
                belowTBL = arr[:, find]
                belowTBL_Tp = T2Tp(belowTBL[0], belowTBL[1])
            tp[zl] = np.mean(belowTBL_Tp[1])
        else:
            LAB1200[zl] = np.nan
            tp[zl] = float(potential_temperature)
            print(LAB1200[zl], tp[zl], locs[:,zl])
        
        arr_pressure=arr.copy().T
        arr_pressure[:,0] -= basement_arr[zl]
        arr_pressure[:,0] = arr_pressure[:,0]/30
        np.savetxt('geoth.dat', arr_pressure, fmt='%12.5f %12.5f')
        # Edit input file
        f=open('input.dat','w')
        f.write('&fit title=\'all_geotherm\',fileinp=\'geoth.dat\',hucrust=%e,hlcrust=%e,vis=%e\n   depthmx=410.,tp=%.0f,mmod=%i,tucrust=%.0f,tlcrust=%.0f,tem=180.,nofj=0,crustk=%.2f\n  ivs=1, filevs=\'lacdevs.dat\',ifit=1,gsize=4.e-3,igrid=1,LON=%e,LAT=%e &end\n'% (h_upper[zl],h_lower[zl],vis[zl],tp[zl],kmod[zl],t[zl]/2.,t[zl]/2.,crustk[zl],locs[zl,0],locs[zl,1]))
        f.close()
        #Run FITPLOT code to find best-fit LAB, fitted geotherms written in multisegment geoth.out file
        cmd = subprocess.Popen('./fit_geotherm', shell=True, stdout=subprocess.PIPE)
        #Extract LAB and misfit from stdout
        cmd2 = subprocess.Popen(('grep "MECH.B.L.=" | tail -n1 |  awk \'{print $2","$4","$6","$8","$12","$14","$16}\' '), shell=True, stdin=cmd.stdout, stdout=subprocess.PIPE)

        output = cmd2.stdout.readline().decode().split(",")
        LAB[zl]=float(output[2])
        MBL[zl]=float(output[0])
        TBL[zl]=float(output[1])
        baseTBL[zl]=float(output[0])+float(output[1])
        misfit[zl]=float(output[3])
        Tmoho[zl]=float(output[4])
        Sflux[zl]=float(output[5])
        Mflux[zl]=float(output[6])

        #subprocess.Popen('cp geoth.out output_geotherms/geoth_'+repr(locs[zl,0])+'_'+repr(locs[zl,1])+'_'+model+'.dat', \
        #shell=True,stdout=subprocess.PIPE)
        arr_out = np.loadtxt('geoth.out')[:,:2]
        arr_out[:,0] += basement_arr[zl]
        np.savetxt(os.path.join(outfold_latest, f'lon_{lon}_lat_{lat}_geoth_6.txt'), arr_out, fmt='%12.5f %12.5f')

        # Remove pre-existing input geotherm
        subprocess.call('rm geoth.dat', shell=True)
        # Print out LAB, misfit, lat and lon
        print (zl, MBL[zl],LAB[zl],TBL[zl],misfit[zl],Tmoho[zl],Sflux[zl],Mflux[zl],locs[zl,0],locs[zl,1])

        t_end=time.time()
        k += 1

main()
exit()

# Mask out NaN regions
test=[None]
mask=[None]

# Read first slice back in
#print(str("Reading in depth %4.0f km" %(seismic_depths[0]/1.e3)))
test = lay_seis_model(model_path = join(seis_pth,\
           str("%ikm.grd" %(seismic_depths[0]/1.e3))));
# Write data
mask=np.reshape(test.data,np.size(test.data))
LAB[mask!=mask]=np.nan
MBL[mask!=mask]=np.nan
TBL[mask!=mask]=np.nan
baseTBL[mask!=mask]=np.nan
misfit[mask!=mask]=np.nan
Tmoho[mask!=mask]=np.nan
Sflux[mask!=mask]=np.nan
Mflux[mask!=mask]=np.nan

HFout=np.zeros((np.shape(temp_arr)[2], 14))
for i in range(np.shape(temp_arr)[2]):
    HFout[i,0]=locs[i,0]
    HFout[i,1]=locs[i,1]
    HFout[i,2]=h_upper[i]*1e+6
    HFout[i,3]=h_lower[i]*1e+6
    HFout[i,4]=crustk[i]
    HFout[i,5]=misfit[i]
    HFout[i,6]=LAB[i]
    HFout[i,7]=MBL[i]
    HFout[i,8]=TBL[i]
    HFout[i,9]=baseTBL[i]
    HFout[i,10]=Tmoho[i]
    HFout[i,11]=Sflux[i]
    HFout[i,12]=Mflux[i]
    HFout[i,13]=tp[i]

format_string=''
for fmt_idx in range(0,4):
    format_string += '%3.1f '
for fmt_idx in range(4,13):
    format_string += '%12.5f'
    format_string += ' '
format_string += '%12.5f'

if os.path.exists('input.dat'):
    subprocess.call('rm input.dat', shell=True)
if os.path.exists('geoth.dat'):
    subprocess.call('rm geoth.dat', shell=True)
if os.path.exists('geoth.out'):
    subprocess.call('rm geoth.out', shell=True)