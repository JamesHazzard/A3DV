import os
import configparser
from datetime import datetime

config_obj = configparser.ConfigParser()
config_obj.read('config.ini')
variables = config_obj["variables"]
folds = config_obj["directories"]

now_date = datetime.strptime(variables["date"], "%Y-%m-%d_%H-%M-%S")
now_date_string = f'RUN_{now_date.strftime("%Y-%m-%d_%H-%M-%S")}'
print(now_date)

potential_temperature = variables["potential_temperature"]
solidus_50km = variables["solidus_50km"]

fold_data_output = os.path.join(folds["data_output"], f'Tp_{potential_temperature}_sol50_{solidus_50km}', 'archive', now_date_string)
fold_data_output_log = os.path.join(folds["data_output"], f'Tp_{potential_temperature}_sol50_{solidus_50km}', 'archive', 'logs')

os.makedirs(fold_data_output, exist_ok=True)
os.makedirs(fold_data_output_log, exist_ok=True)