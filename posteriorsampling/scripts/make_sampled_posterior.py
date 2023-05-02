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
fold_samples_output = os.path.join(fold_data_output, folds["date"], f'Tp_{potential_temperature}_sol50_{solidus_50_km}')

os.makedirs(fold_samples_output, exist_ok=True)

###############################################################################################
# 2. LOAD DATA
###############################################################################################

def load_data():

    data = np.loadtxt(os.path.join(fold_BANCAL22_data, 'samples_postburnin.csv'), skiprows=1)

    return data

data = load_data()

###############################################################################################
# 3. GENERATE RANDOM SAMPLE OF POSTERIOR ANELASTICITY MODELS AND FIND MAP
###############################################################################################

def find_MAP():

    idx_MAP = np.argmax(data[:,0])  # find MAP model by sorting data by posterior probability and extracting the maximum
    data_MAP = data.copy()[idx_MAP,np.arange(0,8,1)]  # get MAP model probability and parameters
    outfile = os.path.join(fold_samples_output, 'MAP_model.txt')
    np.savetxt(outfile, data_MAP, delimiter='\t')  # save MAP model

    return data_MAP

def generate_sample():

    sample_size = 1000
    total_samples = np.shape(data)[0]
    idx_sample = np.random.default_rng().integers(low=0, high=total_samples, size=sample_size)
    data_sample = data.copy()[idx_sample,:]
    data_sample = data_sample[:,np.arange(0,8,1)]
    outfile = os.path.join(fold_samples_output, 'data_sample.txt')
    np.savetxt(outfile, data_sample, delimiter='\t')

    return data_sample

find_MAP()
generate_sample()