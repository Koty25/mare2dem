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
!    vcArg  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      complex(kind=4) cA(10)
      real(kind=4) fB(10)
      real(kind=4) fBha0(10)
      real(kind=4) fBha1(10)
      real(kind=4) fBha2(10)
      real(kind=4) fBla1(10)
      real(kind=4) fBla2(10)
      real(kind=4) fBep1(10)
      real(kind=4) fBep2(10)
      real(kind=4) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      complex(kind=4) cA_I(10*inca)
      real(kind=4) fB_I(10*incb)
      real(kind=4) fBha0_I(10*incb)
      real(kind=4) fBha1_I(10*incb)
      real(kind=4) fBha2_I(10*incb)
      real(kind=4) fBla1_I(10*incb)
      real(kind=4) fBla2_I(10*incb)
      real(kind=4) fBep1_I(10*incb)
      real(kind=4) fBep2_I(10*incb)
      real(kind=4) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      cA( 1)=(-10000.0000,10000.0000)
      cA( 2)=(-7777.7778,7777.7778)
      cA( 3)=(-5555.5557,5555.5557)
      cA( 4)=(-3333.3333,3333.3333)
      cA( 5)=(-1111.1111,1111.1111)
      cA( 6)=(1111.1111,-1111.1111)
      cA( 7)=(3333.3333,-3333.3333)
      cA( 8)=(5555.5557,-5555.5557)
      cA( 9)=(7777.7778,-7777.7778)
      cA(10)=(10000.0000,-10000.0000)
      fB( 1)=2.3561944901923448e+000
      fB( 2)=2.3561944901923448e+000
      fB( 3)=2.3561944901923448e+000
      fB( 4)=2.3561944901923448e+000
      fB( 5)=2.3561944901923448e+000
      fB( 6)=-7.8539816339744828e-001
      fB( 7)=-7.8539816339744828e-001
      fB( 8)=-7.8539816339744828e-001
      fB( 9)=-7.8539816339744828e-001
      fB(10)=-7.8539816339744828e-001

      do i=1,10
          cA_I((i-1)*inca+1)=cA(i)
          fB_I((i-1)*incb+1)=fB(i)
      end do

      call VCARG(vec_len,cA,fBha0)
      call VCARGI(vec_len,cA_I,inca,fBha0_I,incb)

      mode=VML_EP
      call VMCARG(vec_len,cA,fBep1,mode)
      call VMCARGI(vec_len,cA_I,inca,fBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VCARG(vec_len,cA,fBep2)
      call VCARGI(vec_len,cA_I,inca,fBep2_I,incb)

      mode=VML_LA
      call VMCARG(vec_len,cA,fBla1,mode)
      call VMCARGI(vec_len,cA_I,inca,fBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VCARG(vec_len,cA,fBla2)
      call VCARGI(vec_len,cA_I,inca,fBla2_I,incb)

      mode=VML_HA
      call VMCARG(vec_len,cA,fBha1,mode)
      call VMCARGI(vec_len,cA_I,inca,fBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VCARG(vec_len,cA,fBha2)
      call VCARGI(vec_len,cA_I,inca,fBha2_I,incb)

      do i=1,10
          if(fBha0(i) .ne. fBha1(i)) then
              print *,"Error! Difference between VCARG and"
              print *," VMCARG in VML_HA mode detected"
              stop 1
          endif
          if(fBha0_I((i-1)*incb+1) .ne. fBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCARGI and"
              print *," VMCARGI in VML_HA mode detected"
              stop 1
          endif
          if(fBha1(i) .ne. fBha2(i)) then
              print *,"Error! Difference between VCARG and"
              print *," VMCARG in VML_HA mode detected"
              stop 1
          endif
          if(fBha1_I((i-1)*incb+1) .ne. fBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCARGI and"
              print *," VMCARGI in VML_HA mode detected"
              stop 1
          endif
          if(fBla1(i) .ne. fBla2(i)) then
              print *,"Error! Difference between VCARG and"
              print *," VMCARG in VML_LA mode detected"
              stop 1
          endif
          if(fBla1_I((i-1)*incb+1) .ne. fBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCARGI and"
              print *," VMCARGI in VML_LA mode detected"
              stop 1
          endif
          if(fBep1(i) .ne. fBep2(i)) then
              print *,"Error! Difference between VCARG and"
              print *," VMCARG in VML_EP mode detected"
              stop 1
          endif
          if(fBep1_I((i-1)*incb+1) .ne. fBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCARGI and"
              print *," VMCARGI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vcArg test/example program"
      print *,""
      print *,"           Argument                     vcArg"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,cA(i),"    ",fBha0(i)
            CurRMS=srelerr(fB(i),fBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vcArgI test/example program"
      print *,"           Argument                     vcArg"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,cA_I((i-1)*inca+1),"    ",fBha0_I((i-1)*incb+1)
            CurRMS_I=srelerr(fB_I((i-1)*incb+1),fBha0_I((i-1)*incb+1))
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

10    format(E15.7,E15.7,A5,E15.7)
11    format(A,F25.16)

      end
