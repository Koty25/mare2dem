
% This script shows how to take an MT-only data file and a 
% CSEM-only data file and merge them together to create a joint MT-CSEM
% data file using the function m2d_mergeDataFiles.m.
%

sMTDataFile  = 'DemoSynthMT.emdata';
sCSEMDataFile = 'DemoSynthCSEM.emdata';

sOutputFileName = 'DemoSynthCSEM_MT.emdata';

% Output the joint MT-CSEM data file:
m2d_mergeDataFiles(sMTDataFile,sCSEMDataFile,sOutputFileName);
