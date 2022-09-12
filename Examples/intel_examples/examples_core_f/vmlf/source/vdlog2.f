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
!    vdLog2  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=8) :: drelerr

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

      dA(1)=0.1111d0
      dA(2)=1111.2222d0
      dA(3)=2222.3333d0
      dA(4)=3333.4444d0
      dA(5)=4444.5555d0
      dA(6)=5555.6666d0
      dA(7)=6666.7777d0
      dA(8)=7777.8888d0
      dA(9)=8888.9999d0
      dA(10)=9999.0000d0
      dB(1)=-3.1700692781603572d+00
      dB(2)=1.0117931611547215d+01
      dB(3)=1.1117859489416530d+01
      dB(4)=1.1702797948626218d+01
      dB(5)=1.2117823426999088d+01
      dB(6)=1.2439744309294788d+01
      dB(7)=1.2702773906714107d+01
      dB(8)=1.2925162893458976d+01
      dB(9)=1.3117805395452322d+01
      dB(10)=1.3287568102831404d+01

      do i=1,10
          dA_I((i-1)*inca+1)=dA(i)
          dB_I((i-1)*incb+1)=dB(i)
      end do

      call VDLOG2(vec_len,dA,dBha0)
      call VDLOG2I(vec_len,dA_I,inca,dBha0_I,incb)

      mode=VML_EP
      call VMDLOG2(vec_len,dA,dBep1,mode)
      call VMDLOG2I(vec_len,dA_I,inca,dBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VDLOG2(vec_len,dA,dBep2)
      call VDLOG2I(vec_len,dA_I,inca,dBep2_I,incb)

      mode=VML_LA
      call VMDLOG2(vec_len,dA,dBla1,mode)
      call VMDLOG2I(vec_len,dA_I,inca,dBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VDLOG2(vec_len,dA,dBla2)
      call VDLOG2I(vec_len,dA_I,inca,dBla2_I,incb)

      mode=VML_HA
      call VMDLOG2(vec_len,dA,dBha1,mode)
      call VMDLOG2I(vec_len,dA_I,inca,dBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VDLOG2(vec_len,dA,dBha2)
      call VDLOG2I(vec_len,dA_I,inca,dBha2_I,incb)

      do i=1,10
          if(dBha0(i) .ne. dBha1(i)) then
              print *,"Error! Difference between VDLOG2 and"
              print *," VMDLOG2 in VML_HA mode detected"
              stop 1
          endif
          if(dBha0_I((i-1)*incb+1) .ne. dBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDLOG2I and"
              print *," VMDLOG2I in VML_HA mode detected"
              stop 1
          endif
          if(dBha1(i) .ne. dBha2(i)) then
              print *,"Error! Difference between VDLOG2 and"
              print *," VMDLOG2 in VML_HA mode detected"
              stop 1
          endif
          if(dBha1_I((i-1)*incb+1) .ne. dBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDLOG2I and"
              print *," VMDLOG2I in VML_HA mode detected"
              stop 1
          endif
          if(dBla1(i) .ne. dBla2(i)) then
              print *,"Error! Difference between VDLOG2 and"
              print *," VMDLOG2 in VML_LA mode detected"
              stop 1
          endif
          if(dBla1_I((i-1)*incb+1) .ne. dBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDLOG2I and"
              print *," VMDLOG2I in VML_LA mode detected"
              stop 1
          endif
          if(dBep1(i) .ne. dBep2(i)) then
              print *,"Error! Difference between VDLOG2 and"
              print *," VMDLOG2 in VML_EP mode detected"
              stop 1
          endif
          if(dBep1_I((i-1)*incb+1) .ne. dBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VDLOG2I and"
              print *," VMDLOG2I in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vdLog2 test/example program"
      print *,""
      print *,"           Argument                     vdLog2"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA(i),dBha0(i)
            CurRMS=drelerr(dB(i),dBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vdLog2I test/example program"
      print *,"           Argument                     vdLog2"
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
