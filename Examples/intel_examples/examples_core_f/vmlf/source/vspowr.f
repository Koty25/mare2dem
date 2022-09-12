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
!    vsPowr  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: frelerr

      real(kind=4) fA1(10),fA2(10)
      real(kind=4) fB(10)
      real(kind=4) fBha0(10)
      real(kind=4) fBha1(10)
      real(kind=4) fBha2(10)
      real(kind=4) fBla1(10)
      real(kind=4) fBla2(10)
      real(kind=4) fBep1(10)
      real(kind=4) fBep2(10)

      real(kind=8) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      integer, parameter :: incc=7
      real(kind=4) fA1_I(10*inca),fA2_I(10*incb)
      real(kind=4) fB_I(10*incc)
      real(kind=4) fBha0_I(10*incc)
      real(kind=4) fBha1_I(10*incc)
      real(kind=4) fBha2_I(10*incc)
      real(kind=4) fBla1_I(10*incc)
      real(kind=4) fBla2_I(10*incc)
      real(kind=4) fBep1_I(10*incc)
      real(kind=4) fBep2_I(10*incc)
      real(kind=8) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      fA1(1)=0.1000
      fA1(2)=0.8666
      fA1(3)=1.6333
      fA1(4)=2.4000
      fA1(5)=3.1666
      fA1(6)=3.9333
      fA1(7)=4.7000
      fA1(8)=5.4666
      fA1(9)=6.2333
      fA1(10)=7.0000
      fA2(1)=-10.0000;
      fA2(2)=-7.7777;
      fA2(3)=-5.5555;
      fA2(4)=-3.3333;
      fA2(5)=-1.1111;
      fA2(6)=1.1111;
      fA2(7)=3.3333;
      fA2(8)=5.5555;
      fA2(9)=7.7777;
      fA2(10)=10.0000;
      fB(1)=9.9999989760000000e+09
      fB(2)=3.0452833175659180e+00
      fB(3)=6.5510682761669159e-02
      fB(4)=5.4030846804380417e-02
      fB(5)=2.7783754467964172e-01
      fB(6)=4.5796737670898438e+00
      fB(7)=1.7390168762207031e+02
      fB(8)=1.2542625000000000e+04
      fB(9)=1.5173246250000000e+06
      fB(10)=2.8247526400000000e+08

      do i=1,10
          fA1_I((i-1)*inca+1)=fA1(i)
          fA2_I((i-1)*incb+1)=fA2(i)
          fB_I((i-1)*incc+1)=fB(i)
      end do

      call VSPOWR(vec_len,fA1,fA2,fBha0)
      call VSPOWRI(vec_len,fA1_I,inca,fA2_I,incb,fBha0_I,incc)

      mode=VML_EP
      call VMSPOWR(vec_len,fA1,fA2,fBep1,mode)
      call VMSPOWRI(vec_len,fA1_I,inca,fA2_I,incb,fBep1_I,incc,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VSPOWR(vec_len,fA1,fA2,fBep2)
      call VSPOWRI(vec_len,fA1_I,inca,fA2_I,incb,fBep2_I,incc)

      mode=VML_LA
      call VMSPOWR(vec_len,fA1,fA2,fBla1,mode)
      call VMSPOWRI(vec_len,fA1_I,inca,fA2_I,incb,fBla1_I,incc,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VSPOWR(vec_len,fA1,fA2,fBla2)
      call VSPOWRI(vec_len,fA1_I,inca,fA2_I,incb,fBla2_I,incc)

      mode=VML_HA
      call VMSPOWR(vec_len,fA1,fA2,fBha1,mode)
      call VMSPOWRI(vec_len,fA1_I,inca,fA2_I,incb,fBha1_I,incc,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VSPOWR(vec_len,fA1,fA2,fBha2)
      call VSPOWRI(vec_len,fA1_I,inca,fA2_I,incb,fBha2_I,incc)

      do i=1,10
          if(fBha0(i) .ne. fBha1(i)) then
              print *,"Error! Difference between VSPOWR and"
              print *," VMSPOWR in VML_HA mode detected"
              stop 1
          endif
          if(fBha0_I((i-1)*incc+1) .ne. fBha1_I((i-1)*incc+1)) then
              print *,"Error! Difference between VSPOWRI and"
              print *," VMSPOWRI in VML_HA mode detected"
              stop 1
          endif
          if(fBha1(i) .ne. fBha2(i)) then
              print *,"Error! Difference between VSPOWR and"
              print *," VMSPOWR in VML_HA mode detected"
              stop 1
          endif
          if(fBha1_I((i-1)*incc+1) .ne. fBha2_I((i-1)*incc+1)) then
              print *,"Error! Difference between VSPOWRI and"
              print *," VMSPOWRI in VML_HA mode detected"
              stop 1
          endif
          if(fBla1(i) .ne. fBla2(i)) then
              print *,"Error! Difference between VSPOWR and"
              print *," VMSPOWR in VML_LA mode detected"
              stop 1
          endif
          if(fBla1_I((i-1)*incc+1) .ne. fBla2_I((i-1)*incc+1)) then
              print *,"Error! Difference between VSPOWRI and"
              print *," VMSPOWRI in VML_LA mode detected"
              stop 1
          endif
          if(fBep1(i) .ne. fBep2(i)) then
              print *,"Error! Difference between VSPOWR and"
              print *," VMSPOWR in VML_EP mode detected"
              stop 1
          endif
          if(fBep1_I((i-1)*incc+1) .ne. fBep2_I((i-1)*incc+1)) then
              print *,"Error! Difference between VSPOWRI and"
              print *," VMSPOWRI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vsPowr test/example program"
      print *,""
      print *,"           Argument                     vsPowr"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,fA1(i),fA2(i),fBha0(i)
            CurRMS=srelerr(fB(i),fBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vsPowrI test/example program"
      print *,"           Argument                     vsPowr"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,fA1_I((i-1)*inca+1),fA2_I((i-1)*incb+1),           &
     &          fBha0_I((i-1)*incc+1)
            CurRMS_I=srelerr(fB_I((i-1)*incc+1),fBha0_I((i-1)*incc+1))
            if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I
      end do

      print *,""
      if(MaxRMS>=1e-5) then
            print 11,"Error! Relative accuracy is ",MaxRMS
            stop 1
      else
            print 11,"Relative accuracy is ",MaxRMS
      endif

      if(MaxRMS_I>=1e-5) then
            print 11,"Error! Relative strided accuracy is ",MaxRMS_I
            stop 1
      else
            print 11,"Relative strided accuracy is ",MaxRMS_I
      endif

10    format(E25.14,E25.14,E25.14)
11    format(A,F25.16)

      end
