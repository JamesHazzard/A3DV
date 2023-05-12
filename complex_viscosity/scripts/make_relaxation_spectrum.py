import numpy as np
import scipy as sp
import configparser
import os
import time
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from make_loading_history import *
import sys, getopt

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
fold_base = folds["base"]
fold_data_input = folds["data_input"]
fold_data_input_params = os.path.join(fold_data_input, folds["input_anelastic_params"], date, f'Tp_{potential_temperature}_sol50_{solidus_50km}')
fold_data_input_temperature = os.path.join(fold_data_input, folds["input_temperature"], date, f'Tp_{potential_temperature}_sol50_{solidus_50km}')
fold_data_input_viscosity = os.path.join(fold_data_input, folds["input_viscosity"], date, f'Tp_{potential_temperature}_sol50_{solidus_50km}')
fold_data_output = os.path.join(folds["data_output"], date, f'Tp_{potential_temperature}_sol50_{solidus_50km}')
fold_data_output_strain = os.path.join(fold_data_output, 'strain_history')
fold_data_output_viscosity = os.path.join(fold_data_output, 'apparent_viscosity')

os.makedirs(fold_data_output, exist_ok=True)
os.makedirs(fold_data_output_strain, exist_ok=True)
os.makedirs(fold_data_output_viscosity, exist_ok=True)

###############################################################################################
# 2. GET SPATIALLY AVERAGED TEMPERATURE AND VISCOSITY
###############################################################################################

#depth_arr = np.arange(200, 275, 25)
temp_cube = []
visc_cube = []

# initiate variables
fitfile_type = ''
fitfile_index = -1

# read -v variable for type and -i variable for index (if using distribution)
opts, args = getopt.getopt(sys.argv[1:], "v:i:l:s:t:b:")
for opt,arg in opts:
    if opt == '-v':
        fitfile_type = str(arg)
    if opt == '-i':
        fitfile_index = eval(arg)
    if opt == '-l':
        loc = str(arg)
    if opt == '-s':
        stress_model = str(arg)
    if opt == '-t':
        top = eval(arg)
    if opt == '-b':
        bottom = eval(arg)

depth_step = 25
depth_arr = np.arange(top, bottom + depth_step, depth_step) # depth increasing from [top] km to [bottom] km in steps of 25 km

if fitfile_type.strip().lower() in ['map', 'maximum a posteriori', 'most probable']:

    print('MAP model not implemented yet, exiting..')
    exit()

elif (fitfile_type.strip().lower() in ['distribution']) & (fitfile_index >= 0):

    fold_data_output_strain = os.path.join(fold_data_output_strain, 'distribution', 'individual')
    fold_data_output_viscosity = os.path.join(fold_data_output_viscosity, 'distribution', 'individual')
    os.makedirs(fold_data_output_strain, exist_ok=True)
    os.makedirs(fold_data_output_viscosity, exist_ok=True)

    outfile_strain_full_YT16 = os.path.join(fold_data_output_strain, f'strain_full_YT16_location_{loc}_loading_{stress_model}_depthtop_{top}_depthbottom_{bottom}_idx_{fitfile_index}.txt')
    outfile_strain_full_Maxwell = os.path.join(fold_data_output_strain, f'strain_full_Maxwell_location_{loc}_loading_{stress_model}_depthtop_{top}_depthbottom_{bottom}_idx_{fitfile_index}.txt')
    outfile_strain_data_YT16 = os.path.join(fold_data_output_strain, f'strain_data_YT16_location_{loc}_loading_{stress_model}_depthtop_{top}_depthbottom_{bottom}_idx_{fitfile_index}.txt')
    outfile_strain_data_Maxwell = os.path.join(fold_data_output_strain, f'strain_data_Maxwell_location_{loc}_loading_{stress_model}_depthtop_{top}_depthbottom_{bottom}_idx_{fitfile_index}.txt')
    outfile_viscosity = os.path.join(fold_data_output_viscosity, f'apparent_viscosity_location_{loc}_loading_{stress_model}_depthtop_{top}_depthbottom_{bottom}_idx_{fitfile_index}.txt')

    if os.path.exists(outfile_viscosity):

        print("output already exists, exiting..")
        exit()
    
    for idx_depth in range(len(depth_arr)):

        z = depth_arr[idx_depth]
        f_temperature = os.path.join(fold_data_input_temperature,\
            f'location_{loc}', f'depth_{z}_Tp_{potential_temperature}_sol50_{solidus_50km}_idx_{fitfile_index}_{loc}.xyz')
        f_viscosity = os.path.join(fold_data_input_viscosity,\
            f'location_{loc}', f'depth_{z}_Tp_{potential_temperature}_sol50_{solidus_50km}_idx_{fitfile_index}_{loc}.xyz')
        temp_cube.append(np.loadtxt(f_temperature)[:,2])
        #visc_cube.append(np.loadtxt(f_viscosity)[:,2])
    
    f_params = os.path.join(fold_data_input_params, 'data_sample.txt')

else:

    print('No valid user input: please use flags (-v, distribution), (-i, index), (-l, location) and (-s, stress_model)')
    print('index selects the choice of anelasticity model (0-999)')
    print('location selects the region over which to collect temperatures and viscosities')
    print('stress_model selects the loading history to use')
    print('exiting..')
    exit()

# exit if the requested outputs already exist
if os.path.exists(outfile_strain_full_YT16):
        if os.path.exists(outfile_strain_full_Maxwell):
            if os.path.exists(outfile_viscosity):
                print("requested outputs already exist! exiting..")
                #exit()

depth = np.mean(depth_arr)
temp = np.mean(np.concatenate(temp_cube))
#visc = np.mean(np.concatenate(visc_cube))
anyt = np.loadtxt(f_params)[fitfile_index][1:]

###############################################################################################
# 3. DEFINE FUNCTIONS FOR CALCULATING RELAXATION SPECTRUM
###############################################################################################

# Set Yamauchi & Takei 2016 parameters
# Gas constant
R=8.3145
# Reference pressure
Pr=1.5e9
# Reference temperature
TKr=1473.
T0=273.
# Grain size and reference grain size
d=1.e-3
dr=1.e-3
# Reference density
rhor=3291.
# Thermal expansivity
alphaT=3.59e-5
# Bulk modulus
bmod=115.2
# Raw frequency
freq=0.01
Ab=0.664
alphab=0.38
tauP=6.e-5
Teta=0.94
beta=0.
delphi=0.
gamma=5.
lambdaphi=0.

def relaxation_spectrum_YT16(tau, an, temp, depth):    # calculate relaxation spectrum as a function of internal relaxation timescale tau (depends on model params, temperature and pressure)

    mu0=an[0]
    dmudT=an[1]
    dmudP=an[2]
    eta0=10**an[3]
    E=an[4]
    Va=an[5]
    solgrad=an[6]
    sol50=solidus_50km

    TK=temp+273.
    dep=depth
    Pg=(dep/30.)
    P=Pg*1.e9
    Tsol=sol50+(solgrad*(dep-50.))
    Tn=TK/(Tsol+273.)
    
    # Initialise parameters for raw Vs
    if Tn < Teta:
        Aeta=1.
    elif Tn >= Teta and Tn<1.:
        Aeta=np.exp((-1.*((Tn-Teta)/(Tn-(Tn*Teta))))*np.log(gamma))
    else:
        Aeta=(1./gamma)*np.exp(-delphi)
    # Work out viscosity given A
    eta=((eta0*np.exp((E/R)*(1./TK-1./TKr))*np.exp((Va/R)*(P/TK-Pr/TKr)))*Aeta)
    # Unrelaxed compliance
    Ju=1./(1.e9*(mu0+(dmudP*Pg)+(dmudT*(TK-273))))
    tauM=eta*Ju

    if Tn < 0.91:
        Ap=0.01
    elif Tn>=0.91 and Tn<0.96:
        Ap=0.01+(0.4*(Tn-0.91))
    elif Tn>=0.96 and Tn<1.:
        Ap=0.03
    else:
        Ap=0.03+beta
    if Tn < 0.92:
        sigmap=4.
    elif Tn>=0.92 and Tn<1.:
        sigmap=4.+(37.5*(Tn-0.92))
    else:
        sigmap=7.

    taud = tau/tauM
    Xb = Ab*taud**alphab
    Xp = Ap*np.exp(-(((np.log(taud/tauP))**2)/(2*sigmap**2)))
    X = Xb + Xp

    return taud, X, Ju, eta

def relaxation_integrand_YT16(tau, an, temp, depth, t):   # calculate time-dependent integrand required to find creep compliance

    _, X, _, _ = relaxation_spectrum_YT16(tau, an, temp, depth)    # get relaxation spectrum value
    integrand = X*(1 - np.exp(-t/tau))*(1/tau) # multiply relaxation spectrum by exponential term (and divide by tau?) to calculate integrand
    
    return integrand

###############################################################################################
# 4. DEFINE CREEP COMPLIANCE AS A FUNCTION OF RELAXATION SPECTRUM
###############################################################################################

def J_YT16(t, eta, anyt, temp, depth, Ju):

    T_YT16 = np.logspace(start=-256, stop=256, num=3001, base=10.0)
    I_YT16 = cumtrapz(relaxation_integrand_YT16(T_YT16,anyt,temp,depth,t),T_YT16)[-1]
    J_YT16 = Ju + (Ju * I_YT16) + (t / eta) # calculate creep compliance from elastic (instant) + transient (intermediary) + viscous (long timescale) components

    return J_YT16

def J_mx(t, eta, Ju):

    J_mx = Ju + (t / eta) # calculate creep compliance from elastic (instant) + viscous (long timescale) components - Maxwell model contains no transient behaviour

    return J_mx

###############################################################################################
# 5. DEFINE A LOADING HISTORY (RETRIEVE FROM MAKE_LOADING_HISTORY.PY)
###############################################################################################

y2s = 365.25 * 24 * 60 * 60 # year to seconds conversion
models = ["W12", "B18", "Wo15", "I11", "N14", "S20", "S21"]
histories = [load_W12(), load_B18(), load_Wo15(), load_I11(), load_N14(), load_S20(), load_S21()]
mapping = dict(zip(models, histories))
t_history, sigma_history, t_obs_i, t_obs_f, J_obs_elastic = mapping[stress_model]
# remove any duplicate (t, sigma) data points which may occur during construction of loading history
find_dupes = np.where(np.diff(t_history) < 10 / y2s)[0]
t_history = np.delete(t_history.copy(), find_dupes)
sigma_history = np.delete(sigma_history.copy(), find_dupes)
n_timestep = len(t_history)

def plot_loading_history(t, sigma, t_obs_i, t_obs_f, model):

    plt.figure()
    plt.plot(t, sigma, color='r', label='history')
    plt.axvline(t_obs_i, color='k', linestyle='--', label=f'obs start at {t_obs_i}')
    plt.axvline(t_obs_f, color='k', linestyle='--', label=f'obs end at {t_obs_f}')
    plt.ylabel("Stress (arbitrary units)")
    plt.xlabel("Time (years)")
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.savefig(f"testing_{model}_history.jpg", dpi=300)

def plot_all_models():

    models = ["W12", "B18", "Wo15", "I11", "N14", "S20"]
    histories = [load_W12(), load_B18(), load_Wo15(), load_I11(), load_N14(), load_S20()]

    for i in range(len(models)):

        t, sigma, t_i, t_f = histories[i]
        plot_loading_history(t, sigma, t_i, t_f, models[i])

###############################################################################################
# 6. CONVOLVE LOADING HISTORY WITH TIME-DEPENDENT CREEP COMPLIANCE TO CALCULATE STRAIN HISTORY
###############################################################################################

_, _, Ju, eta_ss = relaxation_spectrum_YT16(1e0, anyt, temp, depth)  # extract steady-state parameters Ju and eta_ss
Ju_obs = J_obs_elastic

def calc_transient_strain_YT16():

    ep = np.zeros(n_timestep)
    for i in range(1, n_timestep):
        for j in range(1, i + 1):
            t_conv = y2s * (t_history[i] - t_history[j])    # convolution time is t_conv = t - t'
            J_an = J_YT16(t_conv, eta_ss, anyt, temp, depth, Ju)
            dsig = sigma_history[j] - sigma_history[j-1]
            ep[i] += J_an * dsig
        #print(f"working on time step {i}/{n_timestep-1} of convolution: ep[{i}] is {ep[i]}")

    return ep

def calc_Maxwell_strain(tauM):

    ep_test = np.zeros(n_timestep)
    tauM_test = tauM
    eta_test = (y2s * tauM_test) / Ju
    for i in range(1, n_timestep):
        for j in range(1, i + 1):
            t_conv = y2s * (t_history[i] - t_history[j])
            #J_maxwell = J_mx(t_conv, eta_test, Ju)
            J_maxwell = J_mx(t_conv, eta_test, Ju_obs)
            dsig = sigma_history[j] - sigma_history[j-1]
            ep_test[i] += J_maxwell * dsig

    return eta_test, ep_test

###############################################################################################
# 6. DEFINE RESIDUAL BETWEEN TWO STRAIN HISTORIES
###############################################################################################

def residual_strain(ep_data, ep_test):

    err_ep_data = np.abs(1e-2 * ep_data)    # weight each data point equally by setting a fractional uncertainty of 1%
    chi_squared = np.sum((ep_data - ep_test)**2 / err_ep_data)

    return chi_squared

###############################################################################################
# 7. SIMULATE TRANSIENT AND MAXWELL EVOLUTION (GRID SEARCH OVER TAUM) TO MINIMIZE RESIDUAL
###############################################################################################

ep = calc_transient_strain_YT16()
# strain 'observations' to be compared with Maxwell constitutive law predictions -> interpolate strain profile over obs period
print("interpolating YT16 strain profile..")
spline = sp.interpolate.interp1d(t_history, ep, kind='cubic')
t_fine = np.arange(t_obs_i, t_obs_f, 1 / 365.25)    # assume once daily sampling strategy
t_fine = t_fine[np.where((t_fine >= t_obs_i) & (t_fine <= t_obs_f))]
ep_data = spline(t_fine)
print("done!")

def Maxwell_parameter_sweep(logtauM_lower, logtauM_upper, n_logtauM):

    logtauM = np.linspace(logtauM_lower, logtauM_upper, n_logtauM)
    tauM = 10 ** logtauM
    residual = np.zeros(n_logtauM)

    for idx in range(n_logtauM):

        _, ep_mx = calc_Maxwell_strain(tauM[idx])
        spline_mx = sp.interpolate.interp1d(t_history, ep_mx, kind='cubic')
        ep_mx_data = spline_mx(t_fine)
        residual[idx] = residual_strain(ep_data, ep_mx_data)

    return logtauM, tauM, residual

# coarse tauM sweep
print("running coarse parameter sweep..")
n_coarse = 100
logtauM_coarse_lower = -8
logtauM_coarse_upper = 8
toggle_coarse_sweep = 1

# to make sure we capture the chi-squared minimum, expand the bounds of the coarse parameter sweep if needed!
while toggle_coarse_sweep > 0:

    logtauM_coarse, tauM_coarse, residual_coarse = Maxwell_parameter_sweep(logtauM_coarse_lower, logtauM_coarse_upper, n_coarse)
    idx_minres_coarse = np.argmin(residual_coarse)

    if (idx_minres_coarse > 0) and (idx_minres_coarse < (n_coarse - 1)):
        toggle_coarse_sweep = 0
    elif idx_minres_coarse < 1:
        logtauM_coarse_lower -= (logtauM_coarse_upper - logtauM_coarse_lower)
        logtauM_coarse_upper -= (logtauM_coarse_upper - logtauM_coarse_lower)
    elif idx_minres_coarse > n_coarse - 2:
        logtauM_coarse_lower += (logtauM_coarse_upper - logtauM_coarse_lower)
        logtauM_coarse_upper += (logtauM_coarse_upper - logtauM_coarse_lower)

print("done!")

# plot residual_coarse
#plt.figure()
#plt.plot(logtauM_coarse, np.log10(residual_coarse))
#plt.savefig("test_coarse_residual.jpg", dpi=300)

# fine tauM sweep
n_fine = 100
logtauM_fine_lower = logtauM_coarse[idx_minres_coarse - 1]
logtauM_fine_upper = logtauM_coarse[idx_minres_coarse + 1]
print("running fine parameter sweep..")
logtauM_fine, tauM_fine, residual_fine = Maxwell_parameter_sweep(logtauM_fine_lower, logtauM_fine_upper, n_fine)
print("done!")

# plot residual_fine
#plt.figure()
#plt.plot(logtauM_fine, np.log10(residual_fine))
#plt.savefig("test_fine_residual.jpg", dpi=300)

# pick minimum value
idx_minres_fine = np.argmin(residual_fine)
tauM_optimal = tauM_fine[idx_minres_fine]
eta_mx_optimal, ep_mx_optimal = calc_Maxwell_strain(tauM_optimal)
eta_mx_lower, ep_mx_lower = calc_Maxwell_strain(tauM_fine[idx_minres_fine-1])
eta_mx_upper, ep_mx_upper = calc_Maxwell_strain(tauM_fine[idx_minres_fine+1])
spline_mx_optimal = sp.interpolate.interp1d(t_history, ep_mx_optimal, kind='cubic')
ep_mx_data_optimal = spline_mx_optimal(t_fine)

#print(f"time taken to invert transient viscosity for Maxwell viscosity is {t_end - t_begin} seconds")
#print(f"steady state transient viscosity is {np.log10(eta_ss)}")
#print(f"Maxwell apparent viscosity is {np.log10(eta_mx_lower)}, {np.log10(eta_mx_optimal)}, {np.log10(eta_mx_upper)}")
print(f"(steady-state, Maxwell apparent) viscosity: {np.log10(eta_ss)}, {np.log10(eta_mx_optimal)}")

###############################################################################################
# 8. SAVE CALCULATED STRAIN HISTORIES AND VISCOSITIES TO FILE
###############################################################################################

np.savetxt(outfile_strain_full_YT16, np.stack((t_history, ep)).T, fmt='%12.5f %12e')
np.savetxt(outfile_strain_full_Maxwell, np.stack((t_history, ep_mx_optimal)).T, fmt='%12.5f %12e')
np.savetxt(outfile_strain_data_YT16, np.stack((t_fine, ep_data)).T, fmt='%12.5f %12e')
np.savetxt(outfile_strain_data_Maxwell, np.stack((t_fine, ep_mx_data_optimal)).T, fmt='%12.5f %12e')
np.savetxt(outfile_viscosity, np.array([temp, np.log10(eta_ss), np.log10(eta_mx_optimal)]), fmt = '%12.12f')
t_end = time.time()
print(f"finished in {t_end - t_begin} seconds")