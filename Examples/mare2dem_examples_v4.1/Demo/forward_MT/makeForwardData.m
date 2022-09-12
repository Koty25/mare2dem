%
% This is an example MATLAB code to create a MT data file that
% can be used to generate forward model responses with MARE2DEM.
% 
% Kerry Key
% Lamont-Doherty Earth Observatory
%  

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
% Create a range of frequencies. Here we will use 5 points per decade:
%
nFreqs = 16;
stMT.frequencies = logspace(-4,-1,nFreqs);
 
%
% Set the along profile (y) position of the receivers:
%
nRx     = 25;                       % number of MT receivers 
yRx     = linspace(1000,25000,nRx)';   % this is a column vector
xRx     =  0*ones(nRx,1);           % the x position of MT receivers doesn't matter in 2D

%
% Now we will get the depth of each receiver by interpolating the seafloor
% topography data.
%
% We will subtract 0.1 m from this value to ensure the receiver is slightly
% is located in the ocean not the seafloor.
%
zRx     = interp1(topo(:,1),topo(:,2),yRx)-.1;  
 
%
% To keep this example realistic, we will also set the 2D tilt of the
% receiver (the beta angle) to the seafloor slope:
%
dy         = diff(topo(:,1));
dz         = diff(topo(:,2));
slopeAngle = 180/pi*atan2(dz,dy);
slopePos   = topo(1:end-1,1);  
beta       = interp1(slopePos,slopeAngle,yRx,'previous');



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
% Give the receivers some names:
%
for i = 1:nRx
     stMT.receiverName{i} = sprintf('MT%02i',i);
end

%
% Assemble the MT receiver parameters into the stMT structure:
%
stMT.receivers = [xRx yRx zRx theta alpha beta];

dp = zeros(nFreqs*nRx*4,6);
 

ict = 0;
for ifreq = 1:length(stMT.frequencies)
 
       
    for irx = 1:length(yRx)
         
        dp(ict+1  ,1) = 123;  % log10 Z TE
        dp(ict+1  ,2) = ifreq;
        dp(ict+1  ,3) = irx;
        dp(ict+1  ,4) = irx;

        dp(ict+2  ,1) = 104;  % phase TE
        dp(ict+2  ,2) = ifreq;
        dp(ict+2  ,3) = irx;
        dp(ict+2  ,4) = irx;


        dp(ict+3  ,1) = 125;  % log10 Z TM
        dp(ict+3  ,2) = ifreq;
        dp(ict+3  ,3) = irx;
        dp(ict+3  ,4) = irx;

        dp(ict+4  ,1) = 106;  % phase TM
        dp(ict+4  ,2) = ifreq;
        dp(ict+4  ,3) = irx;
        dp(ict+4  ,4) = irx;

        ict = ict + 4;

  
    end
end
dp = dp(1:ict,:);


%
% Insert into structure to pass to m2d_writeEMData2DFile
%
clear st;
st.stUTM =[];
st.stCSEM = [];
st.comment = 'synth';
st.stMT = stMT;
st.DATA = dp;

%
% Call the data file writing code:
%
m2d_writeEMData2DFile(outputFileName,st) 

%
% That was easy wasn't it?
%
