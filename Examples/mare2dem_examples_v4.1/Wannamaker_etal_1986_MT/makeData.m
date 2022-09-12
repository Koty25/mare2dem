 

topo = load('Hill.txt');

outputFileName = 'Hill.emdata';

stMT = [];
stMT.frequencies = [2 50 2000];
 

% MT Rx

yRx     = [-1500:20:1500]'; 
nRx     = length(yRx);
zRx     = interp1(topo(:,1),topo(:,2),yRx)+.1; % put slightly into the ground
xRx     =  0*ones(nRx,1);
theta   =  0*ones(nRx,1);
alpha   =  0*ones(nRx,1);
beta    =  0*ones(nRx,1);
stMT.receivers = [xRx yRx zRx theta alpha beta];

% Add  receivers with 2D slope parallel E and H sensors
dr = diff(yRx);
dz = diff(zRx);
slope = [0; atan2d(dz,dr)];
stMT.receivers = [stMT.receivers; xRx yRx zRx theta alpha slope];

% Create data array:
dp = zeros(length(stMT.frequencies)*length(yRx)*4*2,6);
 
ict = 0;

% Horizontal receivers:
for ifreq = 1:length(stMT.frequencies)
 
        for irx = 1:length(yRx)
                    
                    itx = irx;
                    
                    dp(ict+1  ,1) = 123;  % log10 Z TE
                    dp(ict+1  ,2) = ifreq;
                    dp(ict+1  ,3) = itx;
                    dp(ict+1  ,4) = irx;

                    dp(ict+2  ,1) = 104;  % phase TE
                    dp(ict+2  ,2) = ifreq;
                    dp(ict+2  ,3) = itx;
                    dp(ict+2  ,4) = irx;

                
                    dp(ict+3  ,1) = 125;  % log10 Z TM
                    dp(ict+3  ,2) = ifreq;
                    dp(ict+3  ,3) = itx;
                    dp(ict+3  ,4) = irx;

                    dp(ict+4  ,1) = 106;  % phase TM
                    dp(ict+4  ,2) = ifreq;
                    dp(ict+4  ,3) = itx;
                    dp(ict+4  ,4) = irx;

                    ict = ict + 4;
                
                    continue
  
    end
end

% Slope parallel receivers:
for ifreq = 1:length(stMT.frequencies)
 
        for irx = 1:length(yRx)
                    
                    itx = irx;
                    
                    dp(ict+1  ,1) = 123;  % log10 Z TE
                    dp(ict+1  ,2) = ifreq;
                    dp(ict+1  ,3) = itx+length(yRx);  % slope parallel H
                    dp(ict+1  ,4) = irx+length(yRx);  % Ex

                    dp(ict+2  ,1) = 104;  % phase TE
                    dp(ict+2  ,2) = ifreq;
                    dp(ict+2  ,3) = itx+length(yRx);  % slope parallel H
                    dp(ict+2  ,4) = irx+length(yRx);  % Ex

                
                    dp(ict+3  ,1) = 125;  % log10 Z TM
                    dp(ict+3  ,2) = ifreq;
                    dp(ict+3  ,3) = itx+length(yRx);  % Hx
                    dp(ict+3  ,4) = irx+length(yRx);  % slope parallel E;

                    dp(ict+4  ,1) = 106;  % phase TM
                    dp(ict+4  ,2) = ifreq;
                    dp(ict+4  ,3) = itx+length(yRx);  % Hx
                    dp(ict+4  ,4) = irx+length(yRx);  % slope parallel E;

                    ict = ict + 4;
                
                    continue
  
    end
end


dp = dp(1:ict,:);

clear st;
st.stUTM = [];
st.stCSEM = [];
st.stMT = stMT;
st.DATA = dp;
m2d_writeEMData2DFile(outputFileName,st) 
