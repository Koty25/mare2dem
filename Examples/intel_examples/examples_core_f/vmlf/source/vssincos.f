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
!    vsSinCos  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      real(kind=4) fA(10)
      real(kind=4) fB(10)
      real(kind=4) fBha0(10)
      real(kind=4) fBha1(10)
      real(kind=4) fBha2(10)
      real(kind=4) fBla1(10)
      real(kind=4) fBla2(10)
      real(kind=4) fBep1(10)
      real(kind=4) fBep2(10)

      real(kind=4) fC(10)
      real(kind=4) fCha0(10)
      real(kind=4) fCha1(10)
      real(kind=4) fCha2(10)
      real(kind=4) fCla1(10)
      real(kind=4) fCla2(10)
      real(kind=4) fCep1(10)
      real(kind=4) fCep2(10)

      real(kind=4) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      integer, parameter :: incc=7
      real(kind=4) fA_I(10*inca)
      real(kind=4) fB_I(10*incb)
      real(kind=4) fBha0_I(10*incb)
      real(kind=4) fBha1_I(10*incb)
      real(kind=4) fBha2_I(10*incb)
      real(kind=4) fBla1_I(10*incb)
      real(kind=4) fBla2_I(10*incb)
      real(kind=4) fBep1_I(10*incb)
      real(kind=4) fBep2_I(10*incb)
      real(kind=4) fC_I(10*incc)
      real(kind=4) fCha0_I(10*incc)
      real(kind=4) fCha1_I(10*incc)
      real(kind=4) fCha2_I(10*incc)
      real(kind=4) fCla1_I(10*incc)
      real(kind=4) fCla2_I(10*incc)
      real(kind=4) fCep1_I(10*incc)
      real(kind=4) fCep2_I(10*incc)
      real(kind=4) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      fA( 1)=-10000.0000
      fA( 2)=-7777.7778
      fA( 3)=-5555.5557
      fA( 4)=-3333.3333
      fA( 5)=-1111.1111
      fA( 6)=1111.1111
      fA( 7)=3333.3333
      fA( 8)=5555.5557
      fA( 9)=7777.7778
      fA(10)=10000.0000
      fB( 1)=3.0561438888825215e-001
      fB( 2)=7.2123131893817349e-001
      fB( 3)=-9.3904866376910323e-001
      fB( 4)=1.0326209316981849e-001
      fB( 5)=8.4827292823844636e-001
      fB( 6)=-8.4827292823844636e-001
      fB( 7)=-1.0326209316981849e-001
      fB( 8)=9.3904866376910323e-001
      fB( 9)=-7.2123131893817349e-001
      fB(10)=-3.0561438888825215e-001
      fC( 1)=-9.5215536825901481e-001
      fC( 2)=6.9269429374198155e-001
      fC( 3)=3.4378424494653864e-001
      fC( 4)=-9.9465418116759940e-001
      fC( 5)=5.2955928772685279e-001
      fC( 6)=5.2955928772685279e-001
      fC( 7)=-9.9465418116759940e-001
      fC( 8)=3.4378424494653864e-001
      fC( 9)=6.9269429374198155e-001
      fC(10)=-9.5215536825901481e-001

      do i=1,10
          fA_I((i-1)*inca+1)=fA(i)
          fB_I((i-1)*incb+1)=fB(i)
          fC_I((i-1)*incc+1)=fC(i)
      end do

      call VSSINCOS(vec_len,fA,fBha0,fCha0)
      call VSSINCOSI(vec_len,fA_I,inca,fBha0_I,incb,fCha0_I,incc)

      mode=VML_EP
      call VMSSINCOS(vec_len,fA,fBep1,fCep1,mode)
      call VMSSINCOSI(vec_len,fA_I,inca,fBep1_I,incb,fCep1_I,incc,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VSSINCOS(vec_len,fA,fBep2,fCep2)
      call VSSINCOSI(vec_len,fA_I,inca,fBep2_I,incb,fCep2_I,incc)

      mode=VML_LA
      call VMSSINCOS(vec_len,fA,fBla1,fCla1,mode)
      call VMSSINCOSI(vec_len,fA_I,inca,fBla1_I,incb,fCla1_I,incc,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VSSINCOS(vec_len,fA,fBla2,fCla2)
      call VSSINCOSI(vec_len,fA_I,inca,fBla2_I,incb,fCla2_I,incc)

      mode=VML_HA
      call VMSSINCOS(vec_len,fA,fBha1,fCha1,mode)
      call VMSSINCOSI(vec_len,fA_I,inca,fBha1_I,incb,fCha1_I,incc,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VSSINCOS(vec_len,fA,fBha2,fCha2)
      call VSSINCOSI(vec_len,fA_I,inca,fBha2_I,incb,fCha2_I,incc)

      do i=1,10
          if((fBha0(i) .ne. fBha1(i)).or.(fCha0(i) .ne. fCha1(i))) then
              print *,"Error! Difference between VSSINCOS and"
              print *," VMSSINCOS in VML_HA mode detected"
              stop 1
          endif
          if((fBha0_I((i-1)*incb+1) .ne. fBha1_I((i-1)*incb+1)) .or.    &
     &        (fCha0_I((i-1)*incc+1) .ne. fCha1_I((i-1)*incc+1))) then
              print *,"Error! Difference between VSSINCOSI and"
              print *," VMSSINCOSI in VML_HA mode detected"
              stop 1
          endif
          if((fBha1(i) .ne. fBha2(i)).or.(fCha1(i) .ne. fCha2(i))) then
              print *,"Error! Difference between VSSINCOS and"
              print *," VMSSINCOS in VML_HA mode detected"
              stop 1
          endif
          if((fBha1_I((i-1)*incb+1) .ne. fBha2_I((i-1)*incb+1)) .or.    &
     &        (fCha1_I((i-1)*incc+1) .ne. fCha2_I((i-1)*incc+1))) then
              print *,"Error! Difference between VSSINCOSI and"
              print *," VMSSINCOSI in VML_HA mode detected"
              stop 1
          endif
          if((fBla1(i) .ne. fBla2(i)).or.(fCla1(i) .ne. fCla2(i))) then
              print *,"Error! Difference between VSSINCOS and"
              print *," VMSSINCOS in VML_LA mode detected"
              stop 1
          endif
          if((fBla1_I((i-1)*incb+1) .ne. fBla2_I((i-1)*incb+1)) .or.    &
     &        (fCla1_I((i-1)*incc+1) .ne. fCla2_I((i-1)*incc+1))) then
              print *,"Error! Difference between VSSINCOSI and"
              print *," VMSSINCOSI in VML_LA mode detected"
              stop 1
          endif
          if((fBep1(i) .ne. fBep2(i)).or.(fCep1(i) .ne. fCep2(i))) then
              print *,"Error! Difference between VSSINCOS and"
              print *," VMSSINCOS in VML_EP mode detected"
              stop 1
          endif
          if((fBep1_I((i-1)*incb+1) .ne. fBep2_I((i-1)*incb+1)) .or.    &
     &        (fCep1_I((i-1)*incc+1) .ne. fCep2_I((i-1)*incc+1))) then
              print *,"Error! Difference between VSSINCOSI and"
              print *," VMSSINCOSI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vsSinCos test/example program"
      print *,""
      print *,"           Argument                           vsSinCos"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,fA(i),fBha0(i),fCha0(i)
            CurRMS=srelerr(fB(i),fBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
            CurRMS=srelerr(fC(i),fCha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vsSinCosI test/example program"
      print *,"           Argument                           vsSinCos"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,fA_I((i-1)*inca+1),fBha0_I((i-1)*incb+1),          &
     &          fCha0_I((i-1)*incc+1)
            CurRMS_I=srelerr(fB_I((i-1)*incb+1),fBha0_I((i-1)*incb+1))
            if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I
            CurRMS_I=srelerr(fC_I((i-1)*incc+1),fCha0_I((i-1)*incc+1))
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
