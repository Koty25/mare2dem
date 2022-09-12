% MATLAB code that reads in a MARE2DEM data file and selectively adds
% static shift(s) to selected stations. This is useful for making synthetic
% data files for testing the "solve static" option of MARE2DEM.
% 
% 
% Run MARE2DEM on the static shifted data file DemoSynthMTstatic.emdata
% add see how well it estimates the static shifts. Then hand edit the data
% file to turn off the solve_static flags for the shifted sites and re-run
% MARE2DEM to see what happens.
%
% Kerry Key
% Scripps Institution of Oceanography
%

%
% The forward response file to read in:
%
fileName = 'DemoSynthMT.emdata';

%
% Specify receivers that will have fake static shifts:
%
iRxStatics = [3 11 19];
nStatic    = log10([15.84 0.166]);  % log10 size of static shifts for TE and TM modes  

%
% Name of the synthetic data file to output:
%
outputFileName = 'DemoSynthMTstatic.emdata';

 
%
% Read the response file in:
% 
st = m2d_readEMData2DFile(fileName);

stMT = st.stMT;
DATA = st.DATA;

%
% Find apparent resistivity data and add the shift:
%
lTEAmp = DATA(:,1) == 123;
lTMAmp = DATA(:,1) == 125;  % data codes for TE and TM log10 ApRes 

for i = 1:length(iRxStatics)
    lRx = DATA(:,4) == iRxStatics(i);
    
    lShiftTE = lTEAmp & lRx;
    lShiftTM = lTMAmp & lRx;
    DATA(lShiftTE,5) = DATA(lShiftTE,5) + nStatic(1);  % since data is log10(apres), we add shift rather than multiply...
    DATA(lShiftTM,5) = DATA(lShiftTM,5) + nStatic(2);
    
end
 
%
% Record static shift flag for these receivers:
%
stMT.receivers(iRxStatics,8) = 1;  

clear stOut;
stOut.stMT = stMT;
stOut.DATA = DATA;
stOut.comment = 'synthetic data with a some static shifts at selected stations';
%
% Write out the synthetic data file:
%
m2d_writeEMData2DFile(outputFileName,stOut) 