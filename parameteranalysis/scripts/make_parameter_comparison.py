import numpy as np

def get_FR20_params():

    anyt = [78.2, -0.020, 2.67, 22.6, 400, 0.092, 0.919]
    err_anyt = [2.2, 0.002, 0.18, 1.6, 288, 5.560, 0.257]

    return anyt, err_anyt

def get_JH23_params():

    anyt = [74.8, -0.0129, 2.04, 23.3, 542, 5.35, 1.63]
    err_anyt = [0.4, 0.00046, 0.06, 0.7, 146, 0.32, 0.14]

    return anyt, err_anyt

def is_overlapping(a1, a2, b1, b2):

    return max(a1, b1) <= min(a2, b2)

x_fr20, y_fr20 = get_FR20_params()
x_jh23, y_jh23 = get_JH23_params()

for i in range(len(x_fr20)):

    if is_overlapping(x_fr20[i] - 2*y_fr20[i], x_fr20[i] + 2*y_fr20[i], x_jh23[i] - 2*y_jh23[i], x_jh23[i] + 2*y_jh23[i]):
        print("True")
    else:
        print(x_fr20[i], y_fr20[i], x_jh23[i], y_jh23[i])