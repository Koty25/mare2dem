
This directory has examples for forward and inverse modeling of MT and CSEM data using MARE2DEM. 

Run MARE2DEM using for example using 5 mpi processes on a laptop:

mpirun -n 5 MARE2DEM Demo.0.resistivity 

When completed, MARE2DEM outputs the response file Demo.0.resp, which can be viewed in plotMARE2DEM_MT.m
or plotMARE2DEM_CSEM.m, depending on the particular data type.

For inversion parameter grids (i.e. models with free parameters) you can add the -S argument 
and MARE2DEM will also output the normalized model parameter sensitivities to files 
Demo.X.sensitivity where X is the inversion iteration number. These can be plotted using 
plotMARE2DEM.m. See the "Overlay" menu.

mpirun -n 5 MARE2DEM -S Demo.0.resistivity 