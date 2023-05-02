import numpy as np
import configparser
import os

###############################################################################################
# 1. GET DIRECTORIES FROM CONFIG
###############################################################################################

config_obj = configparser.ConfigParser()
config_obj.read('config.ini')

folds = config_obj["directories"]
fold_base = folds["base"]
fold_data_input = folds["data_input"]
fold_BANCAL22_runs = os.path.join(fold_data_input, folds["BANCAL22_runs"])
fold_BANCAL22_data = os.path.join(fold_BANCAL22_runs, folds["date"])
potential_temperature = np.loadtxt(os.path.join(fold_BANCAL22_data, 'data', 'potential_temperature', 'potential_temperature.T')).astype(int)
solidus_50_km = np.loadtxt(os.path.join(fold_BANCAL22_data, 'data', 'potential_temperature', 'solidus_50km_temperature.T')).astype(int)
fold_data_output = folds["data_output"]
fold_samples_output = os.path.join(fold_data_output, folds["date"], f'Tp_{potential_temperature}_sol50_{solidus_50_km}', 'parameter_summary')

os.makedirs(fold_samples_output, exist_ok=True)

###############################################################################################
# 2. LOAD DATA
###############################################################################################

def load_data():

    data = np.loadtxt(os.path.join(fold_BANCAL22_data, 'samples_postburnin.csv'), skiprows=1)

    return data

data = load_data()

###############################################################################################
# 3. SUMMARISE DATA
###############################################################################################

p_posterior = data[:,0]
idx_sort = np.argsort(p_posterior)
idx_MAP = idx_sort[-1]
data_MAP = data[idx_MAP, :][0:8]

data_mean = np.nanmean(data, axis=0)[0:8]
data_std = np.nanstd(data, axis=0)[0:8]

np.savetxt(os.path.join(fold_samples_output, 'MAP_model.txt'), data_MAP)
np.savetxt(os.path.join(fold_samples_output, 'mean_model.txt'), data_mean)
np.savetxt(os.path.join(fold_samples_output, 'std_model.txt'), data_std)