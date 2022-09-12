%
% This is an example of MATLAB to code to create a dummy CSEM data file that
% can be used to generate forward model responses with MARE2DEM.
% 
% Kerry Key
% Scripps Institution of Oceanography
%  

minRange = 1000;  % minimum source to receiver distance to model===

%
% Pick a name for the output data file:
%
outputFileName = 'Demo.emdata';

%
% Load in the topography profile since we need it to set the
% receiver depths:
%
topo = load('Topo.prfl');
 
%
% Create a range of frequencies 
% 
stCSEM.frequencies = [0.25 0.75];
 
%
% Set the along profile (y) position of the receivers:
%
nRx     = 25;                        % number of CSEM receivers 
yRx     = linspace(1000,25000,nRx)'; % this is a column vector
xRx     =  0*ones(nRx,1);             

%
% Set the along profile (y) position of the transmitters:
%
nTx     = 25;
yTx     = linspace(0,25000,nTx)';   % this is a column vector
xTx     =  0*ones(nTx,1);                

%
% Now get the depth of each receiver by interpolating the seafloor
% topography data.
%
% We will subtract 0.1 m from this value to ensure the receiver is slightly
% is located in the ocean not the seafloor.  We also use this to set the
% transmitters to be 50 m above the seafloor:
%
zRx     = interp1(topo(:,1),topo(:,2),yRx)-.1; % 0.1 m altitude
zTx     = interp1(topo(:,1),topo(:,2),yTx)-50; % 50 m altitude 

%
% To keep this example realistic, we will also set the 2D tilt of the
% receiver (the beta angle) to the seafloor slope:
%
dy         = diff(topo(:,1));
dz         = diff(topo(:,2));
slopeAngle = 180/pi*atan2(dz,dy);
slopePos   = topo(1:end-1,1);  


if length(slopeAngle) > 1
    beta = interp1(slopePos,slopeAngle,yRx,'previous');
else
    beta = slopeAngle*ones(size(yRx));
end


% Plot Check:
% figure;
% subplot(2,1,1);
% plot(topo(:,1),topo(:,2));axis ij
% ylabel('Depth (m)');
% subplot(2,1,1);
% plot(slopePos,slopeAngle); 
% hold on;
% plot(yRx,beta,'rx')
% ylabel('Slope (degrees clockwise)');

%
% Set the other rotation angles to zero:
%
theta   =  0*ones(nRx,1);
alpha   =  0*ones(nRx,1);


%
% Assemble the CSEM receiver parameters into the stCSEM structure:
%
stCSEM.receivers = [xRx yRx zRx theta alpha beta];

%
% Assemble the CSEM transmitter parameters into the stCSEM structure:
%
stCSEM.transmitters = [xTx yTx zTx 90*ones(size(xTx)) zeros(size(xTx))];  % azimuth (heading) is 90 degrees, so it points along y. Dip = 0;

%
% Give the receivers some names:
%
for i = 1:nRx
    stCSEM.names{i} = sprintf('RX%02i',i);  
end

%
% Give the transmitters some names:
%
for i = 1:nTx
    stCSEM.transmitterName{i} = sprintf('TX%02i',i);  
end


%
% Now apply reciprocity if nTx > nRx
%
if nTx > nRx
    
    A = stCSEM.receivers;
    B = stCSEM.transmitters;
    
    stCSEM.receivers      = B(:,1:3);
    stCSEM.receivers(:,4) = B(:,4)-90; %  Tx azimuth = 90 means Rx x direction = 0, so subtract 90
    stCSEM.receivers(:,5) = 0;         % no alpha angle for transmitters
    stCSEM.receivers(:,6) = B(:,5);    % beta is dip of transmitter
     
    stCSEM.transmitters      = A(:,1:3); 
    stCSEM.transmitters(:,4) = A(:,4)+90; % dipole sources should point along receiver's Ey component
    stCSEM.transmitters(:,5) = A(:,6);    % transmitter dip angle is receiver beta angle
    
    A = stCSEM.transmitterName;
    stCSEM.transmitterName = stCSEM.names;
    stCSEM.names = A;
    
end
    

%
% Set the transmitter type to be electric field dipoles:
%
stCSEM.transmitterType = cell(size(stCSEM.transmitters,1),1);
stCSEM.transmitterType(:) = {'edipole'};


%
% Now create the dummy data array:
%
dp = zeros(length(stCSEM.frequencies)*nRx*nTx*2,6);
 
ict = 0;

for ifreq = 1:length(stCSEM.frequencies)
 
    for itx = 1:size(stCSEM.transmitters,1)
        
    
        for irx = 1:size(stCSEM.receivers,1)
            
            range = abs(stCSEM.receivers(irx,2) - stCSEM.transmitters(itx,2));
            
            if range >= minRange
                
                dp(ict+1  ,1) = 28;  % log10 Amplitude Ey
                dp(ict+1  ,2) = ifreq;
                dp(ict+1  ,3) = itx;
                dp(ict+1  ,4) = irx;

                dp(ict+2  ,1) = 24;  % phase Ey
                dp(ict+2  ,2) = ifreq;
                dp(ict+2  ,3) = itx;
                dp(ict+2  ,4) = irx;
                
                ict = ict + 2;
                
            end        

        end
    end
end
dp = dp(1:ict,:);


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

%
% That was easy wasn't it?
%
