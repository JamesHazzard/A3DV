import numpy as np
from scipy.odr import *
from glob import glob
from os import sys
from os.path import join
import configparser
import os
import matplotlib.pyplot as plt

model = 'ANT-20' # update if using a different seismic model
model_data = 'ant_cont_region_25-400km_downsampled_0.5d'

###############################################################################################
# 1. GET DIRECTORIES FROM CONFIG
###############################################################################################

config_obj = configparser.ConfigParser()
config_obj.read('config.ini')
variables = config_obj["variables"]
potential_temperature = variables["potential_temperature"]
solidus_50km = variables["solidus_50km"]
folds = config_obj["directories"]

fold_latest = os.path.join(folds["data_output"], f'Tp_{potential_temperature}_sol50_{solidus_50km}', 'latest')
fold_LAB1200_latest = os.path.join(fold_latest, folds["output_LAB1200"], 'distribution', 'summary')
fold_baseTBL_latest = os.path.join(fold_latest, folds["output_baseTBL_const_Tp"], 'distribution', 'summary')
fold_LAB1200_to_baseTBL_latest = os.path.join(fold_latest, folds["output_LAB1200_TBL"])

fold_data_input = folds["data_input"]
fold_data_input_LAB_to_baseTBL = os.path.join(fold_data_input, folds["output_LAB1200_TBL"])

os.makedirs(fold_LAB1200_to_baseTBL_latest, exist_ok=True)
os.makedirs(fold_data_input_LAB_to_baseTBL, exist_ok=True)

###############################################################################################
# 2. LOAD LATEST LAB1200 AND BASETBL GRIDS
###############################################################################################

file_LAB1200_latest = os.path.join(fold_LAB1200_latest, 'LAB1200_mean_filtered.xyz')
LAB_lat, LAB_lon, LAB_data = np.loadtxt(file_LAB1200_latest, unpack=True)

file_baseTBL_latest = os.path.join(fold_baseTBL_latest, 'baseTBL_const_Tp_mean_filtered.xyz')
TBL_lat, TBL_lon, TBL_data = np.loadtxt(file_baseTBL_latest, unpack=True)

###############################################################################################
# 3. USE ORTHOGONAL DISTANCE REGRESSION TO FIT RELATIONSHIP BETWEEN LAB1200 AND BASETBL
###############################################################################################

def fitting_func(beta,x):
    return beta[1]*x + beta[0]

data = RealData(LAB_data,TBL_data,50,50)
model = Model(fitting_func)
odr = ODR(data, model, [1,0])
odr.set_job(fit_type=0)
output=odr.run()
output.pprint()
x=LAB_data
y=fitting_func(output.beta,x)

###############################################################################################
# 4. SAVE PLOT OF RELATIONSHIP AND OUTPUT FITTING PARAMETERS
###############################################################################################

plt.figure()
plt.scatter(LAB_data,TBL_data,s=2)
plt.plot(x,y,color='k',linestyle=':',label=r"ODR fit: $y=$"+str(np.round_(output.beta[1],3))+r"$*x+$"+str(np.round_(output.beta[0],3)))
plt.xlabel("LAB1200 depth (km)")
plt.ylabel("baseTBL depth (km)")
plt.legend()
plt.savefig(os.path.join(fold_LAB1200_to_baseTBL_latest, 'scatter.jpg'), dpi=600)

file_output = os.path.join(fold_LAB1200_to_baseTBL_latest, 'LAB1200_to_baseTBL_conversion.txt')
file_output_to_input_dir = os.path.join(fold_data_input_LAB_to_baseTBL, f'LAB1200_to_baseTBL_conversion_Tp_{potential_temperature}_sol50_{solidus_50km}.txt')
output = np.array([output.beta[1],output.beta[0]])
np.savetxt(file_output, output)
np.savetxt(file_output_to_input_dir, output)