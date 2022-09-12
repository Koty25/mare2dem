!===============================================================================
! Copyright 2001-2020 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content:
!    vdPow2o3  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      real(kind=8) dA(10)
      real(kind=8) dB(10)
      real(kind=8) dBha0(10)
      real(kind=8) dBha1(10)
      real(kind=8) dBha2(10)
      real(kind=8) dBla1(10)
      real(kind=8) dBla2(10)
      real(kind=8) dBep1(10)
      real(kind=8) dBep2(10)

      real(kind=8) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      real(kind=8) dA_I(10*inca)
      real(kind=8) dB_I(10*incb)
      real(kind=8) dBha0_I(10*incb)
      real(kind=8) dBha1_I(10*incb)
      real(kind=8) dBha2_I(10*incb)
      real(kind=8) dBla1_I(10*incb)
      real(kind=8) dBla2_I(10*incb)
      real(kind=8) dBep1_I(10*incb)
      real(kind=8) dBep2_I(10*incb)
      real(kind=8) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      dA( 1)=1.0000d+000
      dA( 2)=1.6666d+000
      dA( 3)=2.3333d+000
      dA( 4)=3.0000d+000
      dA( 5)=3.6666d+000
      dA( 6)=4.3333d+000
      dA( 7)=5.0000d+000
      dA( 8)=5.6666d+000
      dA( 9)=6.3333d+000
      dA(10)=7.0000d+000
      dB( 1)=1.0000000000000000d+000
      dB( 2)=1.4056836226901026d+000
      dB( 3)=1.7591939415405531d+000
      dB( 4)=2.0800838230519041d+000
      dB( 5)=2.3778020076302173d+000
      dB( 6)=2.6579440692151435d+000
      dB( 7)=2.9240177382128660d+000
      dB( 8)=3.1784474694102087d+000
      dB( 9)=3.4231035769933422d+000
      dB(10)=3.6593057100229713d+000

      do i=1,10
          dA_I((i-1)*inca+1)=dA(i)
          dB_I((i-1)*incb+1)=dB(i)
      end do

      call VDPOW2O3(vec_len,dA,dBha0)
      call VDPOW2O3I(vec_len,dA_I,inca,dBha0_I,incb)

      mode=VML_EP
      call VMDPOW2O3(vec_len,dA,dBep1,mode)
      call VMDPOW2O3I(vec_len,dA_I,inca,dBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VDPOW2O3(vec_len,dA,dBep2)
      call VDPOW2O3I(vec_len,dA_I,inca,dBep2_I,incb)

      mode=VML_LA
      call VMDPOW2O3(vec_len,dA,dBla1,mode)
      call VMDPOW2O3I(vec_len,dA_I,inca,dBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VDPOW2O3(vec_len,dA,dBla2)
      call VDPOW2O3I(vec_len,dA_I,inca,dBla2_I,incb)

      mode=VML_HA
      call VMDPOW2O3(vec_len,dA,dBha1,mode)
      call VMDPOW2O3I(vec_len,dA_I,inca,dBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VDPOW2O3(vec_len,dA,dBha2)
      call VDPOW2O3I(vec_len,dA_I,inca,dBha2_I,incb)

      do i=1,10
          if(dBha0(i) .ne. dBha1(i)) then
              print *,"Error! Difference between VDPOW2O3 and"
              print *," VMDPOW2O3 in VML_HA mode detected"
              stop 1
          endif
          if(dBha0_I((i-1)*incb+1) .ne. dBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDPOW2O3I and"
              print *," VMDPOW2O3I in VML_HA mode detected"
              stop 1
          endif
          if(dBha1(i) .ne. dBha2(i)) then
              print *,"Error! Difference between VDPOW2O3 and"
              print *," VMDPOW2O3 in VML_HA mode detected"
              stop 1
          endif
          if(dBha1_I((i-1)*incb+1) .ne. dBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDPOW2O3I and"
              print *," VMDPOW2O3I in VML_HA mode detected"
              stop 1
          endif
          if(dBla1(i) .ne. dBla2(i)) then
              print *,"Error! Difference between VDPOW2O3 and"
              print *," VMDPOW2O3 in VML_LA mode detected"
              stop 1
          endif
          if(dBla1_I((i-1)*incb+1) .ne. dBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDPOW2O3I and"
              print *," VMDPOW2O3I in VML_LA mode detected"
              stop 1
          endif
          if(dBep1(i) .ne. dBep2(i)) then
              print *,"Error! Difference between VDPOW2O3 and"
              print *," VMDPOW2O3 in VML_EP mode detected"
              stop 1
          endif
          if(dBep1_I((i-1)*incb+1) .ne. dBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDPOW2O3I and"
              print *," VMDPOW2O3I in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vdPow2o3 test/example program"
      print *,""
      print *,"           Argument                     vdPow2o3"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA(i),dBha0(i)
            CurRMS=drelerr(dB(i),dBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vdPow2o3I test/example program"
      print *,"           Argument                     vdPow2o3"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA_I((i-1)*inca+1),dBha0_I((i-1)*incb+1)
            CurRMS_I=drelerr(dB_I((i-1)*incb+1),dBha0_I((i-1)*incb+1))
            if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I
      end do

      print *,""
      if(MaxRMS>=1e-14) then
            print 11,"Error! Relative accuracy is ",MaxRMS
            stop 1
      else
            print 11,"Relative accuracy is ",MaxRMS
      endif

      if(MaxRMS_I>=1e-14) then
            print 11,"Error! Relative strided accuracy is ",MaxRMS_I
            stop 1
      else
            print 11,"Relative strided accuracy is ",MaxRMS_I
      endif

10    format(E25.14,E25.14)
11    format(A,F25.16)

      end
