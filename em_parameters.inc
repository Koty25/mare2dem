!
! CSEM1D_parameters.inc
!
! Stores magic numbers for CSEM, DC and MT data types
!

!
! Data Type Index Dictionary
!
! These are the codes used in the data file parameter table.
! 
!
! CSEM Data:
! 
! This format is compatible with 2D and 3D CSEM modeling.
! 
!   CSEM Electric Data, Real/Imag:
    integer, parameter :: indRealEx = 1
    integer, parameter :: indImagEx = 2
    integer, parameter :: indRealEy = 3
    integer, parameter :: indImagEy = 4
    integer, parameter :: indRealEz = 5
    integer, parameter :: indImagEz = 6

    
!   CSEM Magnetic Data, Real/Imag:
    integer, parameter :: indRealBx = 11
    integer, parameter :: indImagBx = 12
    integer, parameter :: indRealBy = 13
    integer, parameter :: indImagBy = 14
    integer, parameter :: indRealBz = 15
    integer, parameter :: indImagBz = 16
    
!   CSEM Electric Data, Amp/Phs:
    integer, parameter :: indAmpEx = 21
    integer, parameter :: indPhsEx = 22
    integer, parameter :: indAmpEy = 23
    integer, parameter :: indPhsEy = 24
    integer, parameter :: indAmpEz = 25
    integer, parameter :: indPhsEz = 26
    
!   CSEM Electric Data, log10 Amp:
    integer, parameter :: indLog10AmpEx = 27
    integer, parameter :: indLog10AmpEy = 28
    integer, parameter :: indLog10AmpEz = 29
    
!   CSEM Magnetic Data, Amp/Phs:
    integer, parameter :: indAmpBx = 31
    integer, parameter :: indPhsBx = 32
    integer, parameter :: indAmpBy = 33
    integer, parameter :: indPhsBy = 34
    integer, parameter :: indAmpBz = 35
    integer, parameter :: indPhsBz = 36
    
!   CSEM Magnetic Data, log10 Amp:
    integer, parameter :: indLog10AmpBx = 37
    integer, parameter :: indLog10AmpBy = 38
    integer, parameter :: indLog10AmpBz = 39


!   CSEM Horizontal Polarization Ellipses:
    integer, parameter :: iPEmax = 41
    integer, parameter :: iPEmin = 42
    
    integer, parameter :: iPBmax = 43
    integer, parameter :: iPBmin = 44
    

!   
! Magnetotelluric Data:  
!
! This format is compatible with 2D and 3D MT modeling codes.
!
!   Apparent Resistivity and Phase:
    integer, parameter :: indRhoZXX  = 101   ! reserved for 3D
    integer, parameter :: indPhsZXX  = 102   ! reserved for 3D
    integer, parameter :: indRhoZXY  = 103   ! TE
    integer, parameter :: indPhsZXY  = 104   ! TE  
    integer, parameter :: indRhoZYX  = 105   ! TM
    integer, parameter :: indPhsZYX  = 106   ! TM
    integer, parameter :: indRhoZYY  = 107   ! reserved for 3D
    integer, parameter :: indPhsZYY  = 108   ! reserved for 3D
   
    ! log10(apparent resistivity)
    integer, parameter :: indlog10RhoZXX  = 121 ! reserved for 3D  
    integer, parameter :: indlog10RhoZXY  = 123 ! TE
    integer, parameter :: indlog10RhoZYX  = 125 ! TM
    integer, parameter :: indlog10RhoZYY  = 127 ! reserved for 3D  

        
!   Impedances: 
    integer, parameter :: indRealZXX  = 111   ! reserved for 3D
    integer, parameter :: indImagZXX  = 112   ! reserved for 3D
    integer, parameter :: indRealZXY  = 113   ! TE
    integer, parameter :: indImagZXY  = 114   ! TE
    integer, parameter :: indRealZYX  = 115   ! TM
    integer, parameter :: indImagZYX  = 116   ! TM
    integer, parameter :: indRealZYY  = 117   ! reserved for 3D
    integer, parameter :: indImagZYY  = 118   ! reserved for 3D 
    
! Magnetic Tippers:
    integer, parameter :: indRealMZY  = 133  
    integer, parameter :: indImagMZY  = 134  
    integer, parameter :: indAmpMZY   = 135  
    integer, parameter :: indPhsMZY   = 136  
    
! Raw MT fields:

! For model study purposes only!! These are the calculated fields at the receiver.
! Fields are relative to source field used for boundary conditions and thus should only be 
! compared within a given mode. So Ex,Hy,Hz for TE mode and Hx,Ey,Ez for TM mode. 
!
    integer, parameter :: indRealExMT  = 151 
    integer, parameter :: indImagExMT  = 152  
    integer, parameter :: indRealEyMT  = 153  
    integer, parameter :: indImagEyMT  = 154  
    integer, parameter :: indRealEzMT  = 155  
    integer, parameter :: indImagEzMT  = 156      
    
    integer, parameter :: indRealHxMT  = 161 
    integer, parameter :: indImagHxMT  = 162  
    integer, parameter :: indRealHyMT  = 163  
    integer, parameter :: indImagHyMT  = 164  
    integer, parameter :: indRealHzMT  = 165  
    integer, parameter :: indImagHzMT  = 166     
    
! To do:  other transfer functions ( marine Jz tipper, interstation E-E, H-H ), etc

!   
! DC Resistivity Data:  
! 
integer, parameter :: indAppRes_DC        = 201  
integer, parameter :: indlog10AppRes_DC   = 202  