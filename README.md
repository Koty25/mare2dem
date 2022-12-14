# mare2dem-omeep
Newest MARE2DEM source files customized for the OMEEP project. The original source was directly provided by Kerry key

To compile, run the command line make INCLUDE=./include/omeep.inc
INTEL_PATH=<intel_path> in the repository root dorectory, where <intel_path> is the path to the "intel" directory (where the compilers are installed, it is often /opt/intel).

To compile with score-p (https://www.vi-hps.org/projects/score-p/) automatic instumentation, run the command line make INCLUDE=./include/omeep.inc PREP=scorep INTEL_PATH=<intel_path> in the repository root dorectory, where <intel_path> is the path to the "intel" directory (where the compilers are installed, it is often /opt/intel). Score-p needs to be installed on the machine. More details on installing score-p can be found in its documentation (http://scorepci.pages.jsc.fz-juelich.de/scorep-pipelines/docs/scorep-6.0/pdf/scorep.pdf)

From the original source code download page (https://mare2dem.ucsd.edu/?page_id=108) there are some data set examples. To perform a test execution, go to the inversion_CSEM directory (from the examples root directory) and run the command mpirun -n <nb_processes> <m2d_path>/MARE2DEM Demo.0.resistivity, where <m2d_path> is the path that leads to the MARE2DEM binary, and <nb_processes> is the number of MPI processes. At least two MPI processes are needed to run the application.  
 
The makefile compilation commands have the scorep instrumentation already set. More information about scorep can be found at the link http://www.vi-hps.org/projects/score-p/
