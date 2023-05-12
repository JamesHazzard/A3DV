import numpy as np
import scipy as sp
import configparser
import os
import time
import matplotlib.pyplot as plt
from make_loading_history import *
import sys, getopt
import glob
from sklearn.neighbors import KernelDensity

t_begin = time.time()

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

fold_data_output = os.path.join(folds["data_output"], date, f'Tp_{potential_temperature}_sol50_{solidus_50km}')
fold_data_input_viscosity = os.path.join(fold_data_output, 'apparent_viscosity', 'distribution', 'individual')
fold_data_output_viscosity = os.path.join(fold_data_output, 'apparent_viscosity', 'distribution', 'summary')

os.makedirs(fold_data_output, exist_ok=True)
os.makedirs(fold_data_output_viscosity, exist_ok=True)

###############################################################################################
# 2. GET STEADY STATE AND APPARENT VISCOSITY AND SUMMARISE RESULTS
###############################################################################################

# read -v variable for type and -i variable for index (if using distribution)
opts, args = getopt.getopt(sys.argv[1:], "l:s:t:b:")
for opt,arg in opts:
    if opt == '-l':
        loc = str(arg)
    if opt == '-s':
        stress_model = str(arg)
    if opt == '-t':
        depthtop = eval(arg)
    if opt == '-b':
        depthbottom = eval(arg)

files = glob.glob(os.path.join(fold_data_input_viscosity, f'*location_{loc}_loading_{stress_model}_depthtop_{depthtop}_depthbottom_{depthbottom}_idx_*.txt'))
no_files = len(files)

eta_ss = np.zeros(no_files)
eta_mx = np.zeros(no_files)

for i in range(no_files):

    data = np.loadtxt(files[i])
    eta_ss[i] = data[1]
    eta_mx[i] = data[2]

eta_range = np.linspace(17, 22, 1000)[:, np.newaxis]
sample_range = np.arange(0, no_files, 1)
b_width = 0.035

model_ss = KernelDensity(kernel='gaussian', bandwidth=b_width)
model_ss.fit(eta_ss[:, np.newaxis])
log_dens_ss = model_ss.score_samples(eta_range)

model_mx = KernelDensity(kernel='gaussian', bandwidth=b_width)
model_mx.fit(eta_mx[:, np.newaxis])
log_dens_mx = model_mx.score_samples(eta_range)

np.savetxt(os.path.join(fold_data_output_viscosity, f'apparent_viscosity_location_{loc}_loading_{stress_model}_depthtop_{depthtop}_depthbottom_{depthbottom}_all.txt'), np.stack((eta_ss, eta_mx)).T, fmt='%12.5f')
np.savetxt(os.path.join(fold_data_output_viscosity, f'apparent_viscosity_location_{loc}_loading_{stress_model}_depthtop_{depthtop}_depthbottom_{depthbottom}_summary.txt'), np.array([[np.mean(eta_ss), np.std(eta_ss)], [np.mean(eta_mx), np.std(eta_mx)]]), fmt='%12.5f')
np.savetxt(os.path.join(fold_data_output_viscosity, f'apparent_viscosity_location_{loc}_loading_{stress_model}_depthtop_{depthtop}_depthbottom_{depthbottom}_kde.txt'), np.stack((eta_range[:,0], np.exp(log_dens_ss), np.exp(log_dens_mx))).T, fmt='%21.10f')

print(np.min(eta_ss), np.percentile(eta_ss, 25), np.percentile(eta_ss, 50), np.percentile(eta_ss, 75), np.max(eta_ss))
print(np.min(eta_mx), np.percentile(eta_mx, 25), np.percentile(eta_mx, 50), np.percentile(eta_mx, 75), np.max(eta_mx))
print(f"mean, std over {no_files} models:")
print(np.mean(eta_ss), np.std(eta_ss))
print(np.mean(eta_mx), np.std(eta_mx))

t_end = time.time()
print(f"finished in {t_end - t_begin} seconds")