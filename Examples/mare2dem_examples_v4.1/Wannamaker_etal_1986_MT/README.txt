This is an example showing how to forward model MT responses along topography.

In the data file 'Hill.emdata', the first 151 receivers are idealized MT receivers with horizontal electric and magnetic fields.
These receivers are followed by 151 'tilted' receivers that are located at the same locations as the first 151, but that have 
non-zero beta angles so that y pointing sensors parallel the topographic slope (so they are in the y-z plane). 

The first group of receivers models:
   Z_TE = Ex/Hy
   Z_TM = Ey/Hx

The second group of receivers models:

 Z_TE = Ex / H_parallel 
 Z_TM = E_parallel / Hx.   
 
Although this particular model is a land MT model, the second suite of receivers would be appropriate for marine MT 
where the sensors land on a possibly tilted seafloor.

makeData.m is a MATLAB script that reads in the topo file Hill.txt and creates the horizontal and tilted receivers and then 
writes the data file. 

Run this file in MARE2DEM forward modeling mode using the -F argument, e.g. 

mpirun -np 8 MARE2DEM -F Hill.0.resistivity 

The results can be plotted using the MATLAB script plotMTprofiles.m

This particular model was considered in a seminal paper on topographic effects in MT responses:

Wannamaker, P.E., Stodt, J.A., Rijo, L., 1986. Two-dimensional topographic responses in
magnetotellurics modeled using finite elements. Geophysics 51, 2131–2144.

Check out how the MARE2DEM results compare with their Figures 4-6.

---------------------------------------
Note about simulating real land MT data
---------------------------------------
MARE2DEM now allows you to specify separate receivers for the electric and magnetic field components of an impedance response.
In the data array, electric fields locations are specified by the Rx # column and magnetic fields are specified by the Tx # column.
This can be useful for handling differential tilts between electric and magnetic field components.

In most land MT surveys, the magnetic field sensors are leveled when buried and hence record the horizontal magnetic field, whereas
the electric field dipoles are laid out along the local topography. Therefore, if you are modeling real land MT data in a region 
with significant topography, you can easily setup your data file to handle the 2D aspect of the topography. This is accomplished 
with a slight modification to the example data file here, where the receivers are first listed as horizontal (i.e. no tilts) and
repeated in a second virtual receiver block that has the tilts (i.e. non-zero beta angle, remember its degrees positive down to 
the right).  The data array then uses the horizontal receivers for the magnetic field components and the tilted receivers for the 
electric field components (actually, you only need to do this for the TM mode since it uses E parallel to slope, whereas the TE 
mode uses Ex which is always horizontal, but for simplicity this is ignored below).

 For example the data block would look something like this:
!  Type  Freq #    Tx #    Rx #           Data         StdErr
...
    123       1       1     152              0              0
    104       1       1     152              0              0
    125       1       1     152              0              0
    106       1       1     152              0              0
...
The first two data lines are for log10(TE apparent resistivity) and TE phase while the second two lines are for the TM mode. 
Here the electric fields are taken from receiver #152 (the slope parallel tilted version of receiver #1) whereas the magnetic 
fields are taken from receiver #1 which is horizontal. This has no effect on the TE mode since Z_TE = Ex / Hy. However, for the
TM mode, the slope parallel electric field specified by Rx # 152 means the code will compute Z_TM = E_parallel / Hx,
as collected in the field. 
 
-------------------------------
Note about Hybrid MT Responses:
-------------------------------

With similar syntax to that shown above, MARE2DEM can be used to model hybrid MT responses created from electric and
magnetic field sensors located at different places. In some land MT surveying, it is common to collect continuous electric
field measurements whereas magnetic fields are only recorded at a few locations. This can be modeled with MARE2DEM by
specifying the particular magnetic field receiver location using the Tx # column.  Its all very simple in principle but easy
to mess up, so check and recheck your data formatting code before you trust any MARE2DEM results.


Kerry Key
Scripps Institution of Oceanography
13 April 2013