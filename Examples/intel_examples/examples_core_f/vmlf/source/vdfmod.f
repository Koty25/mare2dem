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
!    vdFmod  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=8) :: drelerr

      real(kind=8) dA1(10),dA2(10)
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
      integer, parameter :: incc=7
      real(kind=8) dA1_I(10*inca),dA2_I(10*incb)
      real(kind=8) dB_I(10*incc)
      real(kind=8) dBha0_I(10*incc)
      real(kind=8) dBha1_I(10*incc)
      real(kind=8) dBha2_I(10*incc)
      real(kind=8) dBla1_I(10*incc)
      real(kind=8) dBla2_I(10*incc)
      real(kind=8) dBep1_I(10*incc)
      real(kind=8) dBep2_I(10*incc)
      real(kind=8) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      dA1(1)=0.1000d0
      dA1(2)=0.8666d0
      dA1(3)=1.6333d0
      dA1(4)=2.4000d0
      dA1(5)=3.1666d0
      dA1(6)=3.9333d0
      dA1(7)=4.7000d0
      dA1(8)=5.4666d0
      dA1(9)=6.2333d0
      dA1(10)=7.0000d0
      dA2(1)=-10.0000d0;
      dA2(2)=-7.7777d0;
      dA2(3)=-5.5555d0;
      dA2(4)=-3.3333d0;
      dA2(5)=-1.1111d0;
      dA2(6)=1.1111d0;
      dA2(7)=3.3333d0;
      dA2(8)=5.5555d0;
      dA2(9)=7.7777d0;
      dA2(10)=10.0000d0;
      dB(1)=1.0000000000000001d-01
      dB(2)=8.6660000000000004d-01
      dB(3)=1.6333000000000000d+00
      dB(4)=2.3999999999999999d+00
      dB(5)=9.4439999999999991d-01
      dB(6)=6.0000000000000009d-01
      dB(7)=1.3667000000000002d+00
      dB(8)=5.4665999999999997d+00
      dB(9)=6.2332999999999998d+00
      dB(10)=7.0000000000000000d+00

      do i=1,10
          dA1_I((i-1)*inca+1)=dA1(i)
          dA2_I((i-1)*incb+1)=dA2(i)
          dB_I((i-1)*incc+1)=dB(i)
      end do

      call VDFMOD(vec_len,dA1,dA2,dBha0)
      call VDFMODI(vec_len,dA1_I,inca,dA2_I,incb,dBha0_I,incc)

      mode=VML_EP
      call VMDFMOD(vec_len,dA1,dA2,dBep1,mode)
      call VMDFMODI(vec_len,dA1_I,inca,dA2_I,incb,dBep1_I,incc,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VDFMOD(vec_len,dA1,dA2,dBep2)
      call VDFMODI(vec_len,dA1_I,inca,dA2_I,incb,dBep2_I,incc)

      mode=VML_LA
      call VMDFMOD(vec_len,dA1,dA2,dBla1,mode)
      call VMDFMODI(vec_len,dA1_I,inca,dA2_I,incb,dBla1_I,incc,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VDFMOD(vec_len,dA1,dA2,dBla2)
      call VDFMODI(vec_len,dA1_I,inca,dA2_I,incb,dBla2_I,incc)

      mode=VML_HA
      call VMDFMOD(vec_len,dA1,dA2,dBha1,mode)
      call VMDFMODI(vec_len,dA1_I,inca,dA2_I,incb,dBha1_I,incc,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VDFMOD(vec_len,dA1,dA2,dBha2)
      call VDFMODI(vec_len,dA1_I,inca,dA2_I,incb,dBha2_I,incc)

      do i=1,10
          if(dBha0(i) .ne. dBha1(i)) then
              print *,"Error! Difference between VDFMOD and"
              print *," VMDFMOD in VML_HA mode detected"
              stop 1
          endif
          if(dBha0_I((i-1)*incc+1) .ne. dBha1_I((i-1)*incc+1)) then
              print *,"Error! Difference between VDFMODI and"
              print *," VMDFMODI in VML_HA mode detected"
              stop 1
          endif
          if(dBha1(i) .ne. dBha2(i)) then
              print *,"Error! Difference between VDFMOD and"
              print *," VMDFMOD in VML_HA mode detected"
              stop 1
          endif
          if(dBha1_I((i-1)*incc+1) .ne. dBha2_I((i-1)*incc+1)) then
              print *,"Error! Difference between VDFMODI and"
              print *," VMDFMODI in VML_HA mode detected"
              stop 1
          endif
          if(dBla1(i) .ne. dBla2(i)) then
              print *,"Error! Difference between VDFMOD and"
              print *," VMDFMOD in VML_LA mode detected"
              stop 1
          endif
          if(dBla1_I((i-1)*incc+1) .ne. dBla2_I((i-1)*incc+1)) then
              print *,"Error! Difference between VDFMODI and"
              print *," VMDFMODI in VML_LA mode detected"
              stop 1
          endif
          if(dBep1(i) .ne. dBep2(i)) then
              print *,"Error! Difference between VDFMOD and"
              print *," VMDFMOD in VML_EP mode detected"
              stop 1
          endif
          if(dBep1_I((i-1)*incc+1) .ne. dBep2_I((i-1)*incc+1)) then
              print *,"Error! Difference between VDFMODI and"
              print *," VMDFMODI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vdFmod test/example program"
      print *,""
      print *,"           Argument                     vdFmod"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA1(i),dA2(i),dBha0(i)
            CurRMS=drelerr(dB(i),dBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vdFmodI test/example program"
      print *,"           Argument                     vdFmod"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,dA1_I((i-1)*inca+1),dA2_I((i-1)*incb+1),           &
     &          dBha0_I((i-1)*incc+1)
            CurRMS_I=drelerr(dB_I((i-1)*incc+1),dBha0_I((i-1)*incc+1))
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
