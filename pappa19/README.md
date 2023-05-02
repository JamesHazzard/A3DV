This directory calculates comparisons between our LAB depth estimates and those of Pappa et al. (2019, JGR).

data_input:
- LAB depth estimates
- polygon files defining continental footprint of Antarctica and west/east divide (see "westvseast" directory)
- difference between LAB depth estimates

plot_output:
- plot comparison between LAB depth estimates

scripts: 
- make_plot_pappa_comparison.sh calculates difference between input LAB depth models
- summary_calculator.py calcultes summary of differences

References:
- Pappa et al. (2019, JGR), DOI: 10.1029/2019JB017997