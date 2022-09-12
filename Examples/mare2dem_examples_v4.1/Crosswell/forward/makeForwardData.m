 
%
% Pick a name for the output data file:
%
outputFileName = 'simple.emdata';

 
%
% Create a range of frequencies 
% 
stCSEM.frequencies = [400];
 
%
% Set the along profile (y) position of the receivers:
%
zRx     = -5+[50:5:290]';
xRx     = 0*zRx;
yRx     = -100 + xRx;
 
theta   = 0*xRx;
alpha   = 0*xRx;
beta    = 0*xRx;
nRx     = length(xRx); 

zTx     = -5+zRx;
xTx     = 0*zTx;
yTx     = 100 + xTx;

azi     = 0*xTx;
dip     = 90 + 0*xTx;
nTx     = length(xTx); 
 

%
% Assemble the CSEM receiver parameters into the stCSEM structure:
%
stCSEM.receivers = [xRx yRx zRx theta alpha beta];

%
% Assemble the CSEM transmitter parameters into the stCSEM structure:
%
stCSEM.transmitters = [xTx yTx zTx azi dip];  % azimuth (heading) is 90 degrees, so it points along y. Dip = 0;

%
% Give the receivers some names:
%
for i = 1:nRx
    stCSEM.receiverName{i} = sprintf('RX%02i',i);  
end

%
% Give the transmitters some names:
%
for i = 1:nTx
    stCSEM.transmitterName{i} = sprintf('TX%02i',i);  
end
 

%
% Set the transmitter type to be electric field dipoles:
%
stCSEM.transmitterType = cell(size(stCSEM.transmitters,1),1);
stCSEM.transmitterType(:) = {'bdipole'};


%
% Now create the dummy data array:
%
dp = zeros(length(stCSEM.frequencies)*nRx*nTx*2,6);
 
ict = 0;

for ifreq = 1:length(stCSEM.frequencies)
 
    for itx = 1:size(stCSEM.transmitters,1)
        
    
        for irx = 1:size(stCSEM.receivers,1)
            
       
                dp(ict+1  ,1) = 39;  % log10 Amplitude Bz
                dp(ict+1  ,2) = ifreq;
                dp(ict+1  ,3) = itx;
                dp(ict+1  ,4) = irx;

                dp(ict+2  ,1) = 36;  % phase Bz
                dp(ict+2  ,2) = ifreq;
                dp(ict+2  ,3) = itx;
                dp(ict+2  ,4) = irx;
                
                ict = ict + 2;
                
          
        end
    end
end
dp = dp(1:ict,:);
dp(:,end) = 1; % dummy error bars

%
% Initialize the output structure:
%
clear st;
st.stUTM =[];
st.stMT = [];
st.stCSEM = stCSEM;
st.stCSEM.phaseConvention = 'lag';
st.DATA = dp;

%
% Call the data file writing code:
%
m2d_writeEMData2DFile(outputFileName,st) 

 