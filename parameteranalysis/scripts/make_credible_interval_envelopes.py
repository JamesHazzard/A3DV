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
fold_data_output_cidf = os.path.join(fold_data_output, folds["date"], f'Tp_{potential_temperature}_sol50_{solidus_50_km}', 'credible_interval_data_fits')
fold_data_output_envelope = os.path.join(fold_data_output_cidf, 'envelope')

os.makedirs(fold_data_output_envelope, exist_ok=True)
viscosity_file = os.path.join(fold_data_output_cidf, 'viscosity', 'viscosity_out.vz')
if os.path.exists(viscosity_file):
    os.remove(viscosity_file)

###############################################################################################
# 2. LOAD DATA
###############################################################################################

def load_data():

    data_loc = fold_data_output_cidf

    plate_loc = data_loc + '/plate'
    adiabat_loc = data_loc + '/adiabat'
    attenuation_loc = data_loc + '/attenuation'
    viscosity_loc = data_loc + '/viscosity'

    plate_full_dir = os.listdir(plate_loc)
    adiabat_full_dir = os.listdir(adiabat_loc)
    attenuation_full_dir = os.listdir(attenuation_loc)
    viscosity_full_dir = os.listdir(viscosity_loc)

    no_models = np.min(np.array([len(plate_full_dir), len(adiabat_full_dir), len(attenuation_full_dir), len(viscosity_full_dir)]))
    print(f'{no_models} to load')

    # load the first model of each data set to find out how many data points reside in each one
    temp_plate = np.loadtxt(plate_loc + '/' + plate_full_dir[0])
    temp_adiabat = np.loadtxt(adiabat_loc + '/' + adiabat_full_dir[0])
    temp_attenuation = np.loadtxt(attenuation_loc + '/' + attenuation_full_dir[0])

    n_plate = np.shape(temp_plate)[0]
    n_adiabat = np.shape(temp_adiabat)[0]   # generalise analogously to plate
    n_attenuation = np.shape(temp_attenuation)[0]  # generalise
    print(n_plate, n_adiabat, n_attenuation)
    
    plate = np.zeros((3, n_plate, no_models))
    adiabat = np.zeros((2, n_adiabat, no_models))   # generalise to n_adiabat
    attenuation = np.zeros((2, n_attenuation, no_models))  # generalise to n_attenuation
    viscosity = np.zeros((2, no_models))    # this is not generalised to number of viscosity data points because inversion not set up for this functionality anyway

    for i in range(no_models):

        print(f"loading data set {i} of {no_models}")

        temp_plate = np.loadtxt(plate_loc + '/' + plate_full_dir[i])
        temp_adiabat = np.loadtxt(adiabat_loc + '/' + adiabat_full_dir[i])
        temp_attenuation = np.loadtxt(attenuation_loc + '/' + attenuation_full_dir[i])
        temp_viscosity = np.loadtxt(viscosity_loc + '/' + viscosity_full_dir[i])

        plate[:,:,i] = temp_plate.T
        adiabat[:,:,i] = temp_adiabat.T
        attenuation[:,:,i] = temp_attenuation.T
        viscosity[:,i] = temp_viscosity.T

    return plate, adiabat, attenuation, viscosity, n_plate, n_adiabat, n_attenuation   

plate, adiabat, attenuation, viscosity, n_plate, n_adiabat, n_attenuation = load_data()

###############################################################################################
# 3. CALCULATE ENVELOPE OF EACH DATA POINT
###############################################################################################

def calculate_distribution_envelope(credible_interval):

    credible_interval/=100
    x_lower = (1-credible_interval)/2
    x_upper = (1+credible_interval)/2

    plate_envelope = np.zeros((5, n_plate))
    adiabat_envelope = np.zeros((4, n_adiabat))
    attenuation_envelope = np.zeros((4, n_attenuation))
    viscosity_envelope = np.zeros((4))

    for i in range(n_plate):

        plate_envelope[0,i] = plate[0,i,0]  # record temperature (fixed) of ith data point
        plate_envelope[1,i] = np.median(plate[1,i,:])   # record median Vs based on distribution of all models
        plate_envelope[2,i] = np.quantile(plate[1,i,:], x_lower)    # x_lower quantile Vs
        plate_envelope[3,i] = np.quantile(plate[1,i,:], x_upper)    # x_upper quantile Vs
        plate_envelope[4,i] = plate[2,i,0]  # record depth assctd. with this data point

    for i in range(n_adiabat):

        adiabat_envelope[3,i] = adiabat[1,i,0]
        adiabat_envelope[0,i] = np.median(adiabat[0,i,:])
        adiabat_envelope[1,i] = np.quantile(adiabat[0,i,:],x_lower)
        adiabat_envelope[2,i] = np.quantile(adiabat[0,i,:],x_upper)

    for i in range(n_attenuation):

        attenuation_envelope[3,i] = attenuation[1,i,0]
        attenuation_envelope[0,i] = np.median(attenuation[0,i,:])
        attenuation_envelope[1,i] = np.quantile(attenuation[0,i,:],x_lower)
        attenuation_envelope[2,i] = np.quantile(attenuation[0,i,:],x_upper)

    viscosity_envelope[3] = viscosity[1,0]
    viscosity_envelope[0] = np.median(viscosity[0,:])
    viscosity_envelope[1] = np.quantile(viscosity[0,:],x_lower)
    viscosity_envelope[2] = np.quantile(viscosity[0,:],x_upper)

    print("envelope calculated")

    return [plate_envelope, adiabat_envelope, attenuation_envelope, viscosity_envelope]

###############################################################################################
# 4. SAVE ENVELOPES
###############################################################################################

def save_envelopes(credible_interval):

    envelope_data = calculate_distribution_envelope(credible_interval)
    save_dir = fold_data_output_envelope

    plate_data = envelope_data[0]
    plate_depths = np.unique(plate_data[4,:])
    no_plate_depths = len(plate_depths)

    for idx_depth in range(no_plate_depths):

        depth = plate_depths[idx_depth]
        temp_data = plate_data[:4,np.where(plate_data[4,:] == depth)[0]]
        
        depth_min = int(depth - 12.5)
        depth_max = int(depth + 12.5)
        
        f_name = os.path.join(save_dir, f'CI_{credible_interval}_envelope_plate_{depth_min}_{depth_max}.TVs')
        np.savetxt(f_name, temp_data.T)
        print(f'plate data for {depth_min}-{depth_max} km layer saved')

    np.savetxt(save_dir+'/CI_'+str(credible_interval)+'_envelope_adiabat.Vsz', envelope_data[1].T)
    print ("adiabat saved")
    np.savetxt(save_dir+'/CI_'+str(credible_interval)+'_envelope_attenuation.Qz', envelope_data[2].T)
    print ("attenuation saved")
    np.savetxt(save_dir+'/CI_'+str(credible_interval)+'_envelope_viscosity.vz', envelope_data[3].T)
    print("viscosity saved")

save_envelopes(99)
save_envelopes(50)
save_envelopes(25)