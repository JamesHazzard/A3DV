all  : fit_geotherm.f
	gfortran fit_geotherm.f -g -O3 -fno-align-commons -fdefault-real-8 -fno-backslash -o fit_geotherm
clean:
	rm fit_geotherm lacdevs.dat output.dat geoth.out
