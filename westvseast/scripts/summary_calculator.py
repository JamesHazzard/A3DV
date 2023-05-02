import numpy as np
import sys, getopt
import scipy.stats as stats
import configparser
import os

###############################################################################################
# 1. LOAD FILE PROVIDED AT INPUT, SUMMARISE, AND SAVE SUMMARY STATISTICS
###############################################################################################

# Load in relevant parameters
opts, args = opts,args=getopt.getopt(sys.argv[1:],"a:")
for opt,arg in opts:
  if opt == '-a':
    root=str(arg)

data = np.loadtxt(root, usecols=(2))
data = data[~np.isnan(data)]
if "eta" in root:
    data = data[np.where(data<=23.0)[0]]
np.savetxt(root+'.summary', np.asarray([np.min(data), np.max(data), np.median(data), stats.median_abs_deviation(data), np.mean(data), np.std(data)]))
