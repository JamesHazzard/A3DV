This diretory sieves the posterior distribution of anelasticity parameters for use in further analysis. 

Data input is taken from BANCAL22 run output directory. Use the config file to adjust where you wish to pull inversion output data from. 

data_output:
- stores *maximum a posteriori (MAP)* set of parameters, and sieved set of parameters

scripts:
- make_sampled_posterior.py script analyses posterior output to find *MAP* and a sieved set of parameters

See Section 4 (including the most relevant Figure 6) of the main manuscript for details of parameter sieving.