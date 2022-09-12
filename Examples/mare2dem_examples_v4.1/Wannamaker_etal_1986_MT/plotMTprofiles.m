
fileName = 'Hill.0.resp';
nRx = 151;

% Read responses into structure st:
st = m2d_readEMData2DFile(fileName);

% Peel off the MT and DATA substructures:
stMT = st.stMT;
DATA =  st.DATA;

iHoriz = DATA(:,3) <= nRx;
iTilt = ~iHoriz;

iTEr = DATA(:,1) == 123;
iTEp = DATA(:,1) == 104;
iTMr = DATA(:,1) == 125;
iTMp = DATA(:,1) == 106;

cc = {'b' 'r' 'g'};
figure;
for ifreq = 1:length(stMT.frequencies)
    
    iuse = DATA(:,2) == ifreq;

    TEr = DATA(iuse&iTEr&iHoriz,7);
    TEp = DATA(iuse&iTEp&iHoriz,7); 
    TMr = DATA(iuse&iTMr&iHoriz,7);
    TMp = DATA(iuse&iTMp&iHoriz,7);     

    iRx = DATA(iuse&iTMp&iHoriz,4);

    TErt = DATA(iuse&iTEr&iTilt,7);
    TEpt = DATA(iuse&iTEp&iTilt,7); 
    TMrt = DATA(iuse&iTMr&iTilt,7);
    TMpt = DATA(iuse&iTMp&iTilt,7);     
    iRxt = DATA(iuse&iTMp&iTilt,4);
    
    subplot(2,2,1);
    hl1(ifreq) = plot( stMT.receivers(iRx,2),TEr,'-','color',cc{ifreq});
    hold on;
    plot(stMT.receivers(iRxt,2),TErt,'--','color',cc{ifreq});
    title ('TE = Ex/Hy (solid = H horizontal, dashed = H parallel )')
    set(gca,'ylim',[1 3])
    
    
    subplot(2,2,3);
    plot( stMT.receivers(iRx,2),TEp,'-','color',cc{ifreq});
    hold on;
    plot( stMT.receivers(iRxt,2),TEpt,'--','color',cc{ifreq});
    set(gca,'ylim',[0 90])
    
    subplot(2,2,2);
    plot( stMT.receivers(iRx,2),TMr,'-','color',cc{ifreq});
    hold on;
    plot( stMT.receivers(iRxt,2),TMrt,'--','color',cc{ifreq});
    title ('TM = Ey/Hx  (solid = E horizontal, dashed = E parallel)')
    set(gca,'ylim',[1 3])
    
    subplot(2,2,4);
    plot( stMT.receivers(iRx,2),TMp,'-','color',cc{ifreq});
    hold on;
    plot(stMT.receivers(iRxt,2),TMpt,'--','color',cc{ifreq});   
    set(gca,'ylim',[0 90])
end
set(gcf,'name','MT responses over a trapezoidal hill, from Wannamaker et al 1986')

subplot(2,2,1);  
legend(hl1,num2str(stMT.frequencies))
ylabel('log10(ohm-m)')

subplot(2,2,3);
xlabel('Position (m)')
ylabel('degrees')

subplot(2,2,4);
xlabel('Position (m)')
ylabel('degrees')

