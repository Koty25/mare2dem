
This folder has all the files already setup to compute the forward response of model Demo.0.resistivity. 

For example using:

mpirun -n 5 MARE2DEM Demo.0.resistivity 

When completed, MARE2DEM outputs the response file Demo.0.resp.

makeForwardData.m is an example script showing how to create the Demo.emdata file that lists the MT 
frequencies, receivers and requested data components.

makeSynthInversionData is an example script showing how to turn the ideal MT responses in Demo.0.resp
into synthetic noisy data that can be used to study how well the MT data and regularized inversion
can recover the true forward model. makeSynthInversionData.m will read in 
the forward response, add random Gaussian noise to the forward response and output the 
synthetic noisy data file to the file DemoSynthMT.emdata, which you can then use to 
in a synthetic inversion run to see how well the data can (and in many places can't) recover the original model.
