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
!    vdAcosh  Example Program Text
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

      dA( 1)=3.5000d+000
      dA( 2)=1114.2222d+000
      dA( 3)=2224.9444d+000
      dA( 4)=3335.6666d+000
      dA( 5)=4446.3888d+000
      dA( 6)=5557.1111d+000
      dA( 7)=6667.8333d+000
      dA( 8)=7778.5555d+000
      dA( 9)=8889.2777d+000
      dA(10)=10000.0000d+000
      dB( 1)=1.9248473002384139d+000
      dB( 2)=7.7090588411869394d+000
      dB( 3)=8.4006343355828150d+000
      dB( 4)=8.8055749765277476d+000
      dB( 5)=9.0929947080081401d+000
      dB( 6)=9.3159808383464178d+000
      dB( 7)=9.4981974184950921d+000
      dB( 8)=9.6522731079815856d+000
      dB( 9)=9.7857482540093379d+000
      dB(10)=9.9034875500361288d+000

      do i=1,10
          dA_I((i-1)*inca+1)=dA(i)
          dB_I((i-1)*incb+1)=dB(i)
      end do

      call VDACOSH(vec_len,dA,dBha0)
      call VDACOSHI(vec_len,dA_I,inca,dBha0_I,incb)

      mode=VML_EP
      call VMDACOSH(vec_len,dA,dBep1,mode)
      call VMDACOSHI(vec_len,dA_I,inca,dBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VDACOSH(vec_len,dA,dBep2)
      call VDACOSHI(vec_len,dA_I,inca,dBep2_I,incb)

      mode=VML_LA
      call VMDACOSH(vec_len,dA,dBla1,mode)
      call VMDACOSHI(vec_len,dA_I,inca,dBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VDACOSH(vec_len,dA,dBla2)
      call VDACOSHI(vec_len,dA_I,inca,dBla2_I,incb)

      mode=VML_HA
      call VMDACOSH(vec_len,dA,dBha1,mode)
      call VMDACOSHI(vec_len,dA_I,inca,dBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VDACOSH(vec_len,dA,dBha2)
      call VDACOSHI(vec_len,dA_I,inca,dBha2_I,incb)

      do i=1,10
          if(dBha0(i) .ne. dBha1(i)) then
              print *,"Error! Difference between VDACOSH and"
              print *," VMDACOSH in VML_HA mode detected"
              stop 1
          endif
          if(dBha0_I((i-1)*incb+1) .ne. dBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDACOSHI and"
              print *," VMDACOSHI in VML_HA mode detected"
              stop 1
          endif
          if(dBha1(i) .ne. dBha2(i)) then
              print *,"Error! Difference between VDACOSH and"
              print *," VMDACOSH in VML_HA mode detected"
              stop 1
          endif
          if(dBha1_I((i-1)*incb+1) .ne. dBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDACOSHI and"
              print *," VMDACOSHI in VML_HA mode detected"
              stop 1
          endif
          if(dBla1(i) .ne. dBla2(i)) then
              print *,"Error! Difference between VDACOSH and"
              print *," VMDACOSH in VML_LA mode detected"
              stop 1
          endif
          if(dBla1_I((i-1)*incb+1) .ne. dBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDACOSHI and"
              print *," VMDACOSHI in VML_LA mode detected"
              stop 1
          endif
          if(dBep1(i) .ne. dBep2(i)) then
              print *,"Error! Difference between VDACOSH and"
              print *," VMDACOSH in VML_EP mode detected"
              stop 1
          endif
          if(dBep1_I((i-1)*incb+1) .ne. dBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDACOSHI and"
              print *," VMDACOSHI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vdAcosh test/example program"
      print *,""
      print *,"           Argument                     vdAcosh"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA(i),dBha0(i)
            CurRMS=drelerr(dB(i),dBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vdAcoshI test/example program"
      print *,"           Argument                     vdAcosh"
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
