import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
y2s = 365.25 * 24 * 60 * 60 # year to seconds conversion
s2y = 1 / y2s
n_tot_global = 75

def calc_timesteps(n_tot, tau_tot, tau_load_i, tau_load_f):

    return int(((tau_load_f - tau_load_i) / tau_tot) * n_tot)

def calc_linear_stress(tau_load_i, tau_load_f, sigma_load_i, sigma_load_f, n_load):

    t_load = np.linspace(tau_load_i, tau_load_f, n_load)
    sigma_load = sigma_load_i + (((t_load - t_load[0]) / (tau_load_f - tau_load_i)) * (sigma_load_f - sigma_load_i))

    return t_load, sigma_load

def load_arbitrary():

    M_elastic = 72.7e+09
    J_elastic = 1./ M_elastic
    
    sigma_max = 1.e+03  # amplitude of stress signal (apparent Maxwell viscosity is independent of this!)
    sigma_i = 1 # initial normalized stress at t=0
    sigma_f = 0 # final normalized stress at t=t_dur
    tau_load = 10 # duration of loading/unloading
    tau_del = 0   # delay between loading/unloading and observation
    #tau_obs_start = tau_load + tau_del  # time when observations started (if tau_del == 0, could be during loading/unloading, otherwise is just tau_laod + tau_del)
    #tau_obs_end = tau_obs_start + 10    # time when observations ended
    tau_obs_start = 10
    tau_obs_end = 20
    tau_tot = tau_obs_end  # total simulation time
    n_timestep = 200   # total number of time steps
    dt = tau_tot / n_timestep
    n_timestep_load = int((tau_load / tau_tot) * n_timestep)  # number of loading/unloading time steps
    n_timestep_del = n_timestep - n_timestep_load   # number of delay time steps

    t_load = np.linspace(0, tau_load, n_timestep_load)   # array corresponding to discretized loading times
    sigma_load = sigma_max * (sigma_i + (sigma_f - sigma_i) * (t_load / tau_load))  # array corresponding to discretized loading stresses

    t_del = np.linspace(t_load[-1], tau_tot, n_timestep_del + 1)[1:]    # array corresponding to discretized delay times
    sigma_del = sigma_max * sigma_f * np.ones(n_timestep_del)   # array corresponding to discretized delay stresses

    t_history = np.concatenate((t_load, t_del)) # array corresponding to discretized full history of times
    idx_observation = np.where((t_history >= tau_obs_start) & (t_history <= tau_obs_end))   # index of data points where observations can be made
    sigma_history = np.concatenate((sigma_load, sigma_del)) # array corresponding to discretized full history of stresses

    return t_history, sigma_history, tau_obs_start, tau_obs_end, J_elastic

def load_W12():

    # Whitehouse et al. (2012)
    M_elastic = 72.7e+09
    J_elastic = 1./ M_elastic

    sigma_max = 100
    n_tot = n_tot_global
    tau_tot = 20e3

    t_history = []
    sigma_history = []

    # linear unloading between 20ka and 5ka (rate 1)
    sigma_load_i = sigma_max
    sigma_load_f = 0.25 * sigma_load_i
    tau_load_i = 0
    tau_load_f = 15e3
    n_load = calc_timesteps(n_tot, tau_tot, tau_load_i, tau_load_f)
    t_load, sigma_load = calc_linear_stress(tau_load_i, tau_load_f, sigma_load_i, sigma_load_f, n_load)
    t_history.append(t_load)
    sigma_history.append(sigma_load)

    # linear unloading between 5ka and 2ka (rate 2)
    sigma_load_i = sigma_load[-1]
    sigma_load_f = 0
    tau_load_i = t_load[-1]
    tau_load_f = tau_load_i + 3e3
    n_load = calc_timesteps(n_tot, tau_tot, tau_load_i, tau_load_f)
    t_load, sigma_load = calc_linear_stress(tau_load_i, tau_load_f, sigma_load_i, sigma_load_f, n_load)
    t_history.append(t_load)
    sigma_history.append(sigma_load)

    # delay between end of deglaciation at 2ka and observations at present day
    tau_load_i = t_load[-1]
    tau_load_f = tau_load_i + 2e3
    n_delay = int(((tau_load_f - tau_load_i) / tau_tot) * n_tot)
    t_delay = np.linspace(tau_load_i, tau_load_f, n_delay)
    sigma_delay = np.zeros(n_delay)
    t_history.append(t_delay)
    sigma_history.append(sigma_delay)

    t_history = np.concatenate(t_history)
    sigma_history = np.concatenate(sigma_history)

    # present day observations (1995-2011) are essentially instantaneous compared to deglaciation history
    tau_obs_i = t_history[-1] - 16
    tau_obs_f = t_history[-1]

    return t_history, sigma_history, tau_obs_i, tau_obs_f, J_elastic

def load_B18():

    # Barletta et al. (2018)
    M_elastic = 72.7e+09
    J_elastic = 1./ M_elastic

    sigma_max = 100
    n_tot = n_tot_global
    tau_tot = 114
    stress_rate_1 = 32.5
    stress_rate_2 = 130
    delta_tau_1 = 102
    delta_tau_2 = 12
    # 0 = sigma_max - (stress_rate_1 / m)*delta_tau_1 - (stress_rate_2 / m)*delta_tau_2 -> solve for normalization m
    m = (1 / sigma_max) * ((stress_rate_1 * delta_tau_1) + (stress_rate_2 * delta_tau_2))

    t_history = []
    sigma_history = []

    # linear unloading between 1900-2002 (rate 1)
    tau_load_i = 0
    tau_load_f = 102
    sigma_load_i = sigma_max
    sigma_load_f = sigma_max - ((stress_rate_1 / m) * (tau_load_f - tau_load_i))
    n_load = calc_timesteps(n_tot, tau_tot, tau_load_i, tau_load_f)
    t_load, sigma_load = calc_linear_stress(tau_load_i, tau_load_f, sigma_load_i, sigma_load_f, n_load)
    t_history.append(t_load)
    sigma_history.append(sigma_load)

    # linear unloading between 2002-2014 (rate 2)
    tau_load_i = t_load[-1]
    tau_load_f = 114
    sigma_load_i = sigma_load[-1]
    sigma_load_f = 0
    n_load = calc_timesteps(n_tot, tau_tot, tau_load_i, tau_load_f)
    t_load, sigma_load = calc_linear_stress(tau_load_i, tau_load_f, sigma_load_i, sigma_load_f, n_load)
    t_history.append(t_load)
    sigma_history.append(sigma_load)

    t_history = np.concatenate(t_history)
    sigma_history = np.concatenate(sigma_history)

    # observations between 2002-2014
    tau_obs_i = 102
    tau_obs_f = 114

    return t_history, sigma_history, tau_obs_i, tau_obs_f, J_elastic

def load_Wo15():

    # Wolstencroft et al. (2015)
    M_elastic = 72.7e+09
    J_elastic = 1./ M_elastic

    # assume W12 deglaciation history but modify observation times
    t_history, sigma_history, _, _, _ = load_W12()
    tau_obs_i = t_history[-1] - 4
    tau_obs_f = t_history[-1]

    return t_history, sigma_history, tau_obs_i, tau_obs_f, J_elastic

def load_S21():

    # Samrat et al. (2021)
    M_elastic = 72.7e+09
    J_elastic = 1./ M_elastic

    sigma_max = 100
    n_tot = n_tot_global
    tau_tot = 21
    stress_rate_1 = 1
    stress_rate_2 = 0.26
    delta_tau_1 = 13
    delta_tau_2 = 8
    # 0 = sigma_max - (stress_rate_1 / m)*delta_tau_1 - (stress_rate_2 / m)*delta_tau_2 -> solve for normalization m
    m = (1 / sigma_max) * ((stress_rate_1 * delta_tau_1) + (stress_rate_2 * delta_tau_2))

    t_history = []
    sigma_history = []

    # linear unloading between 1999-2012 (rate 1)
    tau_load_i = 0
    tau_load_f = 13
    sigma_load_i = sigma_max
    sigma_load_f = sigma_max - ((stress_rate_1 / m) * (tau_load_f - tau_load_i))
    n_load = calc_timesteps(n_tot, tau_tot, tau_load_i, tau_load_f)
    t_load, sigma_load = calc_linear_stress(tau_load_i, tau_load_f, sigma_load_i, sigma_load_f, n_load)
    t_history.append(t_load)
    sigma_history.append(sigma_load)

    # linear unloading between 2012-2020 (rate 2)
    tau_load_i = t_load[-1]
    tau_load_f = 21
    sigma_load_i = sigma_load[-1]
    sigma_load_f = 0
    n_load = calc_timesteps(n_tot, tau_tot, tau_load_i, tau_load_f)
    t_load, sigma_load = calc_linear_stress(tau_load_i, tau_load_f, sigma_load_i, sigma_load_f, n_load)
    t_history.append(t_load)
    sigma_history.append(sigma_load)

    t_history = np.concatenate(t_history)
    sigma_history = np.concatenate(sigma_history)

    # observations between 1999-2020 (GPS site locations each record data over an individual time span - records are incomplete, see Table S3 of Samrat et al., 2021)
    tau_obs_i = 0.01
    tau_obs_f = 20.99

    return t_history, sigma_history, tau_obs_i, tau_obs_f, J_elastic

def load_I11():

    # Ivins et al. (2011)
    M_elastic = 119.25e+09
    J_elastic = 1./ M_elastic

    sigma_max = 100
    n_tot = n_tot_global

    t_load = np.array([-14698.5, -8969.9, -7010.1, -5000, -2990, -2437.24, -1482.46, -829.196, 1452.095808, 1848.502994, 1929.94012, 1992.215569, 1998.203593, 2000.598802, 2002.994012, 2004.191617, 2010])
    sigma_load = np.array([49567.57, 11702.7, 12027.03, 3756.757, 2459.459, 1162.162, 594.5946, 837.8378, 1060, 1116.363636, 792.7272727, 370.909090, 285.454545, 256.363636, 229.090909, 200, 0])

    t_history = np.linspace(t_load[0], t_load[-1], n_tot)
    sigma_history = np.interp(t_history, t_load, sigma_load)
    sigma_history *= (sigma_max / max(sigma_history))
    
    # observations between 2003-2009
    tau_obs_i = 2003
    tau_obs_f = 2009

    return t_history, sigma_history, tau_obs_i, tau_obs_f, J_elastic