This directory allows for the separation of calculated thermomechanical outputs into West vs. East Antarctic continental regions for comparison and summary statistics. 

data_input:
- geotherm grids contains LAB depth and HF grids of choice
- polygons contains continental Antarctic polygon as well as East/West Antarctic polygons, derived from drainage network divides developed by the Goddard Ice Altimetry Group from ICESat data (Zwally et al., 2012)

data_output:
- regionally divided grids and summaries

plot_output:
- histograms and maps plotted based on regionally divided HF/LAB grids

scripts:
- make_divide_grids.sh uses polygons to select regions
- summary_calculator.py calculates simple summary statistics on regional grids

References:
- Zwally et al. (2012), link: https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems