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
!    vdCbrt  Example Program Text
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

      dA( 1)=0.0000d+000
      dA( 2)=1111.1111d+000
      dA( 3)=2222.2222d+000
      dA( 4)=3333.3333d+000
      dA( 5)=4444.4444d+000
      dA( 6)=5555.5555d+000
      dA( 7)=6666.6666d+000
      dA( 8)=7777.7777d+000
      dA( 9)=8888.8888d+000
      dA(10)=10000.0000d+000
      dB( 1)=0.0000000000000000d+000
      dB( 2)=1.0357441651988058d+001
      dB( 3)=1.3049558760397682d+001
      dB( 4)=1.4938015772063830d+001
      dB( 5)=1.6441413774065087d+001
      dB( 6)=1.7710976094006931d+001
      dB( 7)=1.8820720514884833d+001
      dB( 8)=1.9813073109833521d+001
      dB( 9)=2.0714883303976116d+001
      dB(10)=2.1544346900318835d+001

      do i=1,10
          dA_I((i-1)*inca+1)=dA(i)
          dB_I((i-1)*incb+1)=dB(i)
      end do

      call VDCBRT(vec_len,dA,dBha0)
      call VDCBRTI(vec_len,dA_I,inca,dBha0_I,incb)

      mode=VML_EP
      call VMDCBRT(vec_len,dA,dBep1,mode)
      call VMDCBRTI(vec_len,dA_I,inca,dBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VDCBRT(vec_len,dA,dBep2)
      call VDCBRTI(vec_len,dA_I,inca,dBep2_I,incb)

      mode=VML_LA
      call VMDCBRT(vec_len,dA,dBla1,mode)
      call VMDCBRTI(vec_len,dA_I,inca,dBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VDCBRT(vec_len,dA,dBla2)
      call VDCBRTI(vec_len,dA_I,inca,dBla2_I,incb)

      mode=VML_HA
      call VMDCBRT(vec_len,dA,dBha1,mode)
      call VMDCBRTI(vec_len,dA_I,inca,dBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VDCBRT(vec_len,dA,dBha2)
      call VDCBRTI(vec_len,dA_I,inca,dBha2_I,incb)

      do i=1,10
          if(dBha0(i) .ne. dBha1(i)) then
              print *,"Error! Difference between VDCBRT and"
              print *," VMDCBRT in VML_HA mode detected"
              stop 1
          endif
          if(dBha0_I((i-1)*incb+1) .ne. dBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDCBRTI and"
              print *," VMDCBRTI in VML_HA mode detected"
              stop 1
          endif
          if(dBha1(i) .ne. dBha2(i)) then
              print *,"Error! Difference between VDCBRT and"
              print *," VMDCBRT in VML_HA mode detected"
              stop 1
          endif
          if(dBha1_I((i-1)*incb+1) .ne. dBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDCBRTI and"
              print *," VMDCBRTI in VML_HA mode detected"
              stop 1
          endif
          if(dBla1(i) .ne. dBla2(i)) then
              print *,"Error! Difference between VDCBRT and"
              print *," VMDCBRT in VML_LA mode detected"
              stop 1
          endif
          if(dBla1_I((i-1)*incb+1) .ne. dBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDCBRTI and"
              print *," VMDCBRTI in VML_LA mode detected"
              stop 1
          endif
          if(dBep1(i) .ne. dBep2(i)) then
              print *,"Error! Difference between VDCBRT and"
              print *," VMDCBRT in VML_EP mode detected"
              stop 1
          endif
          if(dBep1_I((i-1)*incb+1) .ne. dBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDCBRTI and"
              print *," VMDCBRTI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vdCbrt test/example program"
      print *,""
      print *,"           Argument                     vdCbrt"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA(i),dBha0(i)
            CurRMS=drelerr(dB(i),dBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vdCbrtI test/example program"
      print *,"           Argument                     vdCbrt"
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
