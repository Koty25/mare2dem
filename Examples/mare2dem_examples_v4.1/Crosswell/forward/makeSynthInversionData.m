%
% This is an example of a MATLAB code that reads in a response file, adds
% Gaussian noise to the model responses and then writes out a new
% synthetic data file that can be used for synthetic inversion studies.
% 
% Kerry Key
% Scripps Institution of Oceanography
% 
% Tips:
%
% First run makeFowardData.m to create a dummy data file. Next
% call MARE2DEM using the -F flag to compute forward responses. Then run
% this routine to create the synthetic data file. After this routine you
% need to use Mamba2D to create an inversion model grid, and then you can
% invert the synthetic data using MARE2DEM.
%


%
% The forward response file to read in:
%
fileName = 'simple.0.resp';

%
% Name of the synthetic data file to output:
%
outputFileName = 'simple_synth.emdata';

%
% Specify the amount of Gaussian noise to add (relative):
%
relError = 0.01;  % this means 4 percent random noise.

%
% Read the response file in:
% 
st = m2d_readEMData2DFile(fileName);

% peel off substructures:
stCSEM = st.stCSEM;
DATA = st.DATA;
 

%
% Get log10 ApRes and Phase data:
%
lAmp = DATA(:,1) == 39;  
lPhs = DATA(:,1) == 36;

amp = DATA(lAmp,7); 
phs = DATA(lPhs,7);

%
% Generate synthetic noise:
%
stdA     = 0.4343*relError;  % log10 relative error  
noise    = randn(size(amp));
iLarge   = abs(noise) > 2;   % limit noise to 2 standard deviations, otherwise RMS 1.0 might not be obtainable
noise(iLarge) = sign(noise(iLarge))*2;
noiseAmp = stdA.*noise;

stdP     = 180/pi*relError;   
noise    = randn(size(phs));
iLarge   = abs(noise) > 2;
noise(iLarge) = sign(noise(iLarge))*2;
noisePhs = stdP.*noise;

%
% Add noise to the model responses:
%
outA = amp + noiseAmp;
outP = phs + noisePhs;
 
%
%  Create a new data array with the synthetic data:
%
dp = DATA(:,1:4);
dp(lAmp,5) = outA;
dp(lPhs,5) = outP;

dp(lAmp,6) = stdA;
dp(lPhs,6) = stdP;

% insert into a strucutre:
clear stOut;
stOut.stCSEM = stCSEM;
stOut.DATA = dp;


%
% Write out the synthetic data file:
%
m2d_writeEMData2DFile(outputFileName,stOut) 