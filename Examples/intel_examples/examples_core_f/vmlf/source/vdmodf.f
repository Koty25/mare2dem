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
!    vdModf  Example Program Text
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

      real(kind=8) dC(10)
      real(kind=8) dCha0(10)
      real(kind=8) dCha1(10)
      real(kind=8) dCha2(10)
      real(kind=8) dCla1(10)
      real(kind=8) dCla2(10)
      real(kind=8) dCep1(10)
      real(kind=8) dCep2(10)

      real(kind=8) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      integer, parameter :: incc=7
      real(kind=8) dA_I(10*inca)
      real(kind=8) dB_I(10*incb)
      real(kind=8) dBha0_I(10*incb)
      real(kind=8) dBha1_I(10*incb)
      real(kind=8) dBha2_I(10*incb)
      real(kind=8) dBla1_I(10*incb)
      real(kind=8) dBla2_I(10*incb)
      real(kind=8) dBep1_I(10*incb)
      real(kind=8) dBep2_I(10*incb)
      real(kind=8) dC_I(10*incc)
      real(kind=8) dCha0_I(10*incc)
      real(kind=8) dCha1_I(10*incc)
      real(kind=8) dCha2_I(10*incc)
      real(kind=8) dCla1_I(10*incc)
      real(kind=8) dCla2_I(10*incc)
      real(kind=8) dCep1_I(10*incc)
      real(kind=8) dCep2_I(10*incc)
      real(kind=8) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      dA( 1)=-2.0000d+000
      dA( 2)=-1.5555d+000
      dA( 3)=-1.1111d+000
      dA( 4)=-0.6666d+000
      dA( 5)=-0.2222d+000
      dA( 6)=0.2222d+000
      dA( 7)=0.6666d+000
      dA( 8)=1.1111d+000
      dA( 9)=1.5555d+000
      dA(10)=2.0000d+000
      dB( 1)=-2.0000000000000000d+000
      dB( 2)=-1.0000000000000000d+000
      dB( 3)=-1.0000000000000000d+000
      dB( 4)=0.0000000000000000d+000
      dB( 5)=0.0000000000000000d+000
      dB( 6)=0.0000000000000000d+000
      dB( 7)=0.0000000000000000d+000
      dB( 8)=1.0000000000000000d+000
      dB( 9)=1.0000000000000000d+000
      dB(10)=2.0000000000000000d+000
      dC( 1)=0.0000000000000000d+000
      dC( 2)=-5.5550000000000010d-001
      dC( 3)=-1.1109999999999998d-001
      dC( 4)=-6.6659999999999997d-001
      dC( 5)=-2.2220000000000001d-001
      dC( 6)=2.2220000000000001d-001
      dC( 7)=6.6659999999999997d-001
      dC( 8)=1.1109999999999998d-001
      dC( 9)=5.5550000000000010d-001
      dC(10)=0.0000000000000000d+000

      do i=1,10
          dA_I((i-1)*inca+1)=dA(i)
          dB_I((i-1)*incb+1)=dB(i)
          dC_I((i-1)*incc+1)=dC(i)
      end do

      call VDMODF(vec_len,dA,dBha0,dCha0)
      call VDMODFI(vec_len,dA_I,inca,dBha0_I,incb,dCha0_I,incc)

      mode=VML_EP
      call VMDMODF(vec_len,dA,dBep1,dCep1,mode)
      call VMDMODFI(vec_len,dA_I,inca,dBep1_I,incb,dCep1_I,incc,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VDMODF(vec_len,dA,dBep2,dCep2)
      call VDMODFI(vec_len,dA_I,inca,dBep2_I,incb,dCep2_I,incc)

      mode=VML_LA
      call VMDMODF(vec_len,dA,dBla1,dCla1,mode)
      call VMDMODFI(vec_len,dA_I,inca,dBla1_I,incb,dCla1_I,incc,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VDMODF(vec_len,dA,dBla2,dCla2)
      call VDMODFI(vec_len,dA_I,inca,dBla2_I,incb,dCla2_I,incc)

      mode=VML_HA
      call VMDMODF(vec_len,dA,dBha1,dCha1,mode)
      call VMDMODFI(vec_len,dA_I,inca,dBha1_I,incb,dCha1_I,incc,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VDMODF(vec_len,dA,dBha2,dCha2)
      call VDMODFI(vec_len,dA_I,inca,dBha2_I,incb,dCha2_I,incc)

      do i=1,10
          if((dBha0(i) .ne. dBha1(i)).or.(dCha0(i) .ne. dCha1(i))) then
              print *,"Error! Difference between VDMODF and"
              print *," VMDMODF in VML_HA mode detected"
              stop 1
          endif
          if((dBha0_I((i-1)*incb+1) .ne. dBha1_I((i-1)*incb+1)) .or.    &
     &        (dCha0_I((i-1)*incc+1) .ne. dCha1_I((i-1)*incc+1))) then
              print *,"Error! Difference between VDMODFI and"
              print *," VMDMODFI in VML_HA mode detected"
              stop 1
          endif
          if((dBha1(i) .ne. dBha2(i)).or.(dCha1(i) .ne. dCha2(i))) then
              print *,"Error! Difference between VDMODF and"
              print *," VMDMODF in VML_HA mode detected"
              stop 1
          endif
          if((dBha1_I((i-1)*incb+1) .ne. dBha2_I((i-1)*incb+1)) .or.    &
     &        (dCha1_I((i-1)*incc+1) .ne. dCha2_I((i-1)*incc+1))) then
              print *,"Error! Difference between VDMODFI and"
              print *," VMDMODFI in VML_HA mode detected"
              stop 1
          endif
          if((dBla1(i) .ne. dBla2(i)).or.(dCla1(i) .ne. dCla2(i))) then
              print *,"Error! Difference between VDMODF and"
              print *," VMDMODF in VML_LA mode detected"
              stop 1
          endif
          if((dBla1_I((i-1)*incb+1) .ne. dBla2_I((i-1)*incb+1)) .or.    &
     &        (dCla1_I((i-1)*incc+1) .ne. dCla2_I((i-1)*incc+1))) then
              print *,"Error! Difference between VDMODFI and"
              print *," VMDMODFI in VML_LA mode detected"
              stop 1
          endif
          if((dBep1(i) .ne. dBep2(i)).or.(dCep1(i) .ne. dCep2(i))) then
              print *,"Error! Difference between VDMODF and"
              print *," VMDMODF in VML_EP mode detected"
              stop 1
          endif
          if((dBep1_I((i-1)*incb+1) .ne. dBep2_I((i-1)*incb+1)) .or.    &
     &        (dCep1_I((i-1)*incc+1) .ne. dCep2_I((i-1)*incc+1))) then
              print *,"Error! Difference between VDMODFI and"
              print *," VMDMODFI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vdModf test/example program"
      print *,""
      print *,"           Argument                             vdModf"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA(i),dBha0(i),dCha0(i)
            CurRMS=drelerr(dB(i),dBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
            CurRMS=drelerr(dC(i),dCha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vdModfI test/example program"
      print *,"           Argument                             vdModf"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA_I((i-1)*inca+1),dBha0_I((i-1)*incb+1),          &
     &          dCha0_I((i-1)*incc+1)
            CurRMS_I=drelerr(dB_I((i-1)*incb+1),dBha0_I((i-1)*incb+1))
            if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I
            CurRMS_I=drelerr(dC_I((i-1)*incc+1),dCha0_I((i-1)*incc+1))
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

10    format(E25.14,E25.14,E25.14)
11    format(A,F25.16)

      end
