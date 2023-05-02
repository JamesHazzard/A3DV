import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os 
import configparser

###############################################################################################
# 1. GET DIRECTORIES FROM CONFIG
###############################################################################################

config_obj = configparser.ConfigParser()
config_obj.read('config.ini')
folds = config_obj["directories"]
fold_base = folds["base"]
fold_data_input = folds["data_input"]
fold_input_magmatism = os.path.join(fold_data_input, folds["input_magmatism"])
fold_input_LAB_models = os.path.join(fold_data_input, folds["input_LAB_models"])
fold_input_conductive_isotherms = os.path.join(fold_data_input, folds["input_conductive_isotherms"])
fold_data_output = folds["data_output"]

###############################################################################################
# 2. READ IN LAB-MAGMATISM JOINT DATA SETS AND USE MONTE CARLO SIM TO EVALUATE CORRELATION
###############################################################################################

def monte_carlo_spearman(model):    # use uncertainties to simulate range of data sets and calculate correlation

    print(f"simulating range of data for model {model}")
    data = np.loadtxt(os.path.join(fold_data_output, f'Antarctica_magmatism_LAB_{model}_blockmean.xyaszs'), usecols=(2,3,4,5))
    data_filter = data[np.where((data[:,1]<10.0))[0],:]
    trials = 100000
    s_all = np.zeros(trials)
    p_all = np.zeros(trials)
    s_filter = np.zeros(trials)
    p_filter = np.zeros(trials)
    for i in range(trials):
        x = data[:,0] + np.random.normal(scale=np.abs(data[:,1]))
        y = data[:,2] + np.random.normal(scale=np.abs(data[:,3]))
        stat = stats.spearmanr(x, y)
        s_all[i] = stat[0]
        p_all[i] = stat[1]
        x = data_filter[:,0] + np.random.normal(scale=np.abs(data_filter[:,1]))
        y = data_filter[:,2] + np.random.normal(scale=np.abs(data_filter[:,3]))
        stat = stats.spearmanr(x, y)
        s_filter[i] = stat[0]
        p_filter[i] = stat[1]
    np.savetxt(os.path.join(fold_data_output, f'Antarctica_magmatism_LAB_{model}_spearman.spsp'), np.stack((s_all, p_all, s_filter, p_filter)).T)

def summarise_spearman(model):  # load and print correlations

    file_path = os.path.join(fold_data_output, f'Antarctica_magmatism_LAB_{model}_spearman.spsp')
    skip_rerun=True
    if (os.path.exists(file_path)) & (skip_rerun):
        s1, p1, s2, p2 = np.loadtxt(file_path, unpack=True)
    else:
        monte_carlo_spearman(model)
        s1, p1, s2, p2 = np.loadtxt(file_path, unpack=True)

    #print("s is", np.mean(s2), np.std(s2), np.mean(s2)-np.std(s2), np.mean(s2)+np.std(s2))
    #print("p is", np.median(p2), stats.median_abs_deviation(p2))
    return s2, p2

def record_spearman():  # save correlations

    s = np.zeros((3, 100000))
    p = np.zeros((3, 100000))
    s[0],p[0]=summarise_spearman('ANT-20')
    s[1],p[1]=summarise_spearman('SL2013')
    s[2],p[2]=summarise_spearman('CAM2016')
    np.savetxt(os.path.join(fold_data_output, 'Antarctica_magmatism_LAB_spearman_coefficient.txt'), np.stack((np.mean(s,axis=1), np.std(s,axis=1))).T, header='ANT-20, SL2013, CAM2016')
    np.savetxt(os.path.join(fold_data_output, 'Antarctica_magmatism_LAB_spearman_p_value.txt'), np.stack((np.median(p,axis=1), stats.median_abs_deviation(p,axis=1))).T, header='ANT-20, SL2013, CAM2016')

def histogram_plot():   # (optional) plot histogram

    s = np.zeros((3, 100000))
    p = np.zeros((3, 100000))
    s[0],p[0]=summarise_spearman('ANT-20')
    s[1],p[1]=summarise_spearman('SL2013')
    s[2],p[2]=summarise_spearman('CAM2016')
    fig,ax=plt.subplots()
    plt.hist(s[0],bins=100, alpha=0.5, density=True, label='ANT-20')
    plt.hist(s[1],bins=100, alpha=0.5, density=True, label='SL2013')
    plt.hist(s[2],bins=100, alpha=0.5, density=True, label='CAM2016')
    #ax.set_xlim(0,0.6)
    #ax.set_xlabel("p-value")
    ax.axvline(0.296, linestyle='--')
    ax.set_xlabel("Spearman's Rank Correlation Coefficient")
    ax.set_ylabel("Probability Density")
    plt.legend(loc="upper left")
    plt.show()

###############################################################################################
# 3. RUN
###############################################################################################

record_spearman()
