import numpy as np
import pandas as pd
import csv
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
fold_input_magmatism = os.path.join(fold_data_input, folds["input_magmatism"])
fold_input_LAB_models = os.path.join(fold_data_input, folds["input_LAB_models"])
fold_input_conductive_isotherms = os.path.join(fold_data_input, folds["input_conductive_isotherms"])
fold_data_output = folds["data_output"]

###############################################################################################
# 2. DEFINE FUNCTIONS TO CALCULATE AGE OF SAMPLE BASED ON COMBINATION OF LABELS/VALUES
###############################################################################################

def calc_age_from_label(age_label):

    era_labels = ['PROTEROZOIC', 'PALEOZOIC', 'MESOZOIC', 'CENOZOIC']
    era_age_max = [2500, 541, 251.902, 66]
    era_age_min = [541, 251.902, 66, 0]
    period_labels = ['CAMBRIAN', 'ORDOVICIAN', 'SILURIAN', 'DEVONIAN', 'CARBONIFEROUS', 'PERMIAN', 'TRIASSIC', 'JURASSIC', 'CRETACEOUS', 'PALEOGENE', 'NEOGENE', 'QUATERNARY']
    period_age_max = [541, 485.4, 443.8, 419.2, 358.9, 298.9, 251.902, 201.3, 145, 66, 23.03, 2.58]
    period_age_min = [485.4, 443.8, 419.2, 358.9, 298.9, 251.902, 201.3, 145, 66, 23.03, 2.58, 0]
    epoch_labels = ['PALEOCENE', 'EOCENE', 'OLIGOCENE', 'MIOCENE', 'PLIOCENE', 'PLEISTOCENE', 'HOLOCENE', 'PLIO-PLEISTOCENE']
    epoch_age_max = [66, 56, 33.9, 23.03, 5.333, 2.58, 0.0117, 5.333]
    epoch_age_min = [56, 33.9, 23.03, 5.333, 2.58, 0.0117, 0, 0.0117]
    deep_time_labels = era_labels + period_labels + epoch_labels
    deep_time_age_min = era_age_min + period_age_min + epoch_age_min
    deep_time_age_max = era_age_max + period_age_max + epoch_age_max

    iter = 0
    k = 0
    age_min = np.NaN
    age_max = np.NaN
    while iter < 0.5:
        if age_label.strip() in deep_time_labels[k]:
            age_min = deep_time_age_min[k]
            age_max = deep_time_age_max[k]
            iter = 1
        k+=1

    return age_min, age_max

def process_label(x):

    epoch_labels = ['PALEOCENE', 'EOCENE', 'OLIGOCENE', 'MIOCENE', 'PLIOCENE', 'PLEISTOCENE', 'HOLOCENE', 'PLIO-PLEISTOCENE']
    period_labels = ['PALEOGENE', 'NEOGENE', 'QUATERNARY']
    era_labels = ['CENOZOIC']
    #input is of the form e.g. ['CENOZOIC [xxx]', 'PLIOCENE [zzz]', 'NEOGENE [yyy]'] (not necessarily in tree order)
    age_epoch = []
    age_period = []
    age_era = []
    try:
        length = len(x)
        for j in range(length):
            x[j] = x[j].split()[0]
            if x[j] in epoch_labels:
                age_epoch.append(calc_age_from_label(x[j]))
            elif x[j] in period_labels:
                age_period.append(calc_age_from_label(x[j]))
            elif x[j] in era_labels:
                age_era.append(calc_age_from_label(x[j]))
        if len(age_epoch) > 0:
            data = np.asarray(age_epoch)
            arg = np.where(data==np.min(data))[0][0]
            age_min = data[arg][0]
            age_max = data[arg][1]
        elif len(age_period) > 0:
            data = np.asarray(age_period)
            arg = np.where(data==np.min(data))[0][0]
            age_min = data[arg][0]
            age_max = data[arg][1]
        elif len(age_era) > 0:
            data = np.asarray(age_era)
            arg = np.where(data==np.min(data))[0][0]
            age_min = data[arg][0]
            age_max = data[arg][1]
        else:
            age_min = np.NaN
            age_max = np.NaN
    except:
        age_min = np.NaN
        age_max = np.NaN

    return age_min, age_max

def process_label_min_age(x):

    age_min, age_max = process_label(x)

    return age_min

def process_label_max_age(x):

    age_min, age_max = process_label(x)

    return age_max

def process_min_age(x):

    age = []
    try:
        length = len(x)
        for j in range(length):
            x[j]=x[j].split()[0]
            age.append(float(x[j]))
        x = min(age)/(10**6)
    except:
        x = np.NaN

    return x

def process_max_age(x):

    age = []
    try:
        length = len(x)
        for j in range(length):
            x[j]=x[j].split()[0]
            age.append(float(x[j]))
        x = max(age)/(10**6)
    except:
        x = np.NaN
        
    return x

###################################################################################################
# 3. DEFINE FUNCTIONS TO EXTRACT DATA FOR EACH SPECIFIC DATA SET (SARBAS, 2008 & BALL ET AL., 2021)
###################################################################################################

def extract_GEOROC_data():

    print("Extracting magmatism data from Sarbas (2008) aka GEOROC..")
    # Load in the GEOROC data base and remove the rows where no age data at all is present
    file_input_path = os.path.join(fold_input_magmatism, "ANTARCTICA_rift.csv")
    file_output_path = os.path.join(fold_data_output, "ANTARCTICA_rift_filtered.csv")
    df = pd.read_csv(file_input_path, delimiter=',',\
        usecols=['LATITUDE MIN', 'LATITUDE MAX', 'LONGITUDE MIN', 'LONGITUDE MAX', 'MIN. AGE (YRS.)', 'MAX. AGE (YRS.)', 'AGE'])\
        .dropna(subset=['MIN. AGE (YRS.)', 'MAX. AGE (YRS.)', 'AGE'], how='all')
    # Split geological time labels for each sample (row) into a list (e.g. "Cenozoic / Oligocene" -> ['Cenozoic', 'Oligocene'])
    df['AGE'] = df['AGE'].str.split(" / ")
    # Find minimum and maximum age of sample based on implied age of geological time label(s) [where present]
    df['MIN. AGE (FROM LABEL, Ma)'] = df['AGE'].map(lambda x: process_label_min_age(x))
    df['MAX. AGE (FROM LABEL, Ma)'] = df['AGE'].map(lambda x: process_label_max_age(x))
    # Find minimum and maximum age of sample based on quantitative age record [where present]
    df['MIN. AGE (Ma)'] = df['MIN. AGE (YRS.)'].str.split("/").map(lambda x: process_min_age(x))
    df['MAX. AGE (Ma)'] = df['MAX. AGE (YRS.)'].str.split("/").map(lambda x: process_max_age(x))
    # Convert age columns into numerical values to allow use of np.where function
    df[['MIN. AGE (Ma)', 'MAX. AGE (Ma)', 'MIN. AGE (FROM LABEL, Ma)', 'MAX. AGE (FROM LABEL, Ma)']] = \
        df[['MIN. AGE (Ma)', 'MAX. AGE (Ma)', 'MIN. AGE (FROM LABEL, Ma)', 'MAX. AGE (FROM LABEL, Ma)']].apply(pd.to_numeric)
    # Estimate age from average of quantitative min/max where present, and from geological time period implied min/max where present [if neither present, set as np.NaN]
    df['NUMERICAL (0) OR LABEL (1) FLAG'] = \
        np.where(df['MIN. AGE (Ma)'].isnull() | df['MAX. AGE (Ma)'].isnull(),\
        np.where(df['MIN. AGE (FROM LABEL, Ma)'].isnull() | df['MAX. AGE (FROM LABEL, Ma)'].isnull(), np.NaN, int(1)), int(0))
    df['AGE ESTIMATE (Ma)'] = np.where(df['MIN. AGE (Ma)'].isnull() | df['MAX. AGE (Ma)'].isnull(),\
        np.where(df['MIN. AGE (FROM LABEL, Ma)'].isnull() | df['MAX. AGE (FROM LABEL, Ma)'].isnull(), np.NaN, (df['MIN. AGE (FROM LABEL, Ma)'] + df['MAX. AGE (FROM LABEL, Ma)']) / 2), (df['MIN. AGE (Ma)'] + df['MAX. AGE (Ma)']) / 2)
    # Estimate age uncertainty from half-range of quantitative min/max where present, and from geological time period implied min/max where present [if neither present, set as np.NaN]
    df['AGE ESTIMATE UNCERTAINTY (Ma)'] = np.where(df['MIN. AGE (Ma)'].isnull() | df['MAX. AGE (Ma)'].isnull(),\
        np.where(df['MIN. AGE (FROM LABEL, Ma)'].isnull() | df['MAX. AGE (FROM LABEL, Ma)'].isnull(), np.NaN, -(df['MIN. AGE (FROM LABEL, Ma)'] - df['MAX. AGE (FROM LABEL, Ma)']) / 2), -(df['MIN. AGE (Ma)'] - df['MAX. AGE (Ma)']) / 2)
    df = df[df['LATITUDE MIN'] >= -85.0]
    df_save = df[['LATITUDE MIN', 'LATITUDE MAX', 'LONGITUDE MIN', 'LONGITUDE MAX', 'AGE ESTIMATE (Ma)', 'AGE ESTIMATE UNCERTAINTY (Ma)', 'NUMERICAL (0) OR LABEL (1) FLAG']]\
        .dropna(subset=['AGE ESTIMATE (Ma)'])
    df_save.to_csv(file_output_path, index=False)
    print("done")

    return None

def extract_BALL_data():

    print("Extracting magmatism data from Ball et al. (2021)..")
    file_input_path = os.path.join(fold_input_magmatism, "BALL_2021_Antarctica.csv")
    file_output_path = os.path.join(fold_data_output, "BALL_filtered.csv")
    df = pd.read_csv(file_input_path, delimiter=',', usecols=['Latitude', 'Longitude', 'Age', 'Uncertainty'])
    df_save = df[['Latitude', 'Longitude', 'Age', 'Uncertainty']].dropna(subset=['Age'])
    df_save.to_csv(file_output_path, index=False)
    print("done")

    return None

###############################################################################################
# 4. RUN PROCESS
###############################################################################################

extract_GEOROC_data()
extract_BALL_data()
