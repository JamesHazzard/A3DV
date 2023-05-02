cool.o  : RHCW18_plate_model.f90 global.f90
	gfortran -o cool global.f90 RHCW18_plate_model.f90 -O3 -ftree-vectorize -g
clean:
	rm -f cool ../../data_output/hf-*.dat ../../data_output/Tgrid-*.dat ../../data_output/depth-*.dat