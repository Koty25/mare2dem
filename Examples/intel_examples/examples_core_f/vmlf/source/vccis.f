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
!    vcCIS  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      real(kind=4) fA(10)
      complex(kind=4) cB(10)
      complex(kind=4) cBha0(10)
      complex(kind=4) cBha1(10)
      complex(kind=4) cBha2(10)
      complex(kind=4) cBla1(10)
      complex(kind=4) cBla2(10)
      complex(kind=4) cBep1(10)
      complex(kind=4) cBep2(10)

      real(kind=4) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      real(kind=4) fA_I(10*inca)
      complex(kind=4) cB_I(10*incb)
      complex(kind=4) cBha0_I(10*incb)
      complex(kind=4) cBha1_I(10*incb)
      complex(kind=4) cBha2_I(10*incb)
      complex(kind=4) cBla1_I(10*incb)
      complex(kind=4) cBla2_I(10*incb)
      complex(kind=4) cBep1_I(10*incb)
      complex(kind=4) cBep2_I(10*incb)
      real(kind=4) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      fA( 1)=-10000.0000
      fA( 2)=-7777.7777
      fA( 3)=-5555.5555
      fA( 4)=-3333.3333
      fA( 5)=-1111.1111
      fA( 6)=1111.1111
      fA( 7)=3333.3333
      fA( 8)=5555.5555
      fA( 9)=7777.7777
      fA(10)=10000.0000
      cB( 1)=(-9.5215536825901481e-001,3.0561438888825215e-001)
      cB( 2)=(6.9269429374198155e-001,7.2123131893817349e-001)
      cB( 3)=(3.4378424494653864e-001,-9.3904866376910323e-001)
      cB( 4)=(-9.9465418116759940e-001,1.0326209316981849e-001)
      cB( 5)=(5.2955928772685279e-001,8.4827292823844636e-001)
      cB( 6)=(5.2955928772685279e-001,-8.4827292823844636e-001)
      cB( 7)=(-9.9465418116759940e-001,-1.0326209316981849e-001)
      cB( 8)=(3.4378424494653864e-001,9.3904866376910323e-001)
      cB( 9)=(6.9269429374198155e-001,-7.2123131893817349e-001)
      cB(10)=(-9.5215536825901481e-001,-3.0561438888825215e-001)

      do i=1,10
          fA_I((i-1)*inca+1)=fA(i)
          cB_I((i-1)*incb+1)=cB(i)
      end do

      call VCCIS(vec_len,fA,cBha0)
      call VCCISI(vec_len,fA_I,inca,cBha0_I,incb)

      mode=VML_EP
      call VMCCIS(vec_len,fA,cBep1,mode)
      call VMCCISI(vec_len,fA_I,inca,cBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VCCIS(vec_len,fA,cBep2)
      call VCCISI(vec_len,fA_I,inca,cBep2_I,incb)

      mode=VML_LA
      call VMCCIS(vec_len,fA,cBla1,mode)
      call VMCCISI(vec_len,fA_I,inca,cBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VCCIS(vec_len,fA,cBla2)
      call VCCISI(vec_len,fA_I,inca,cBla2_I,incb)

      mode=VML_HA
      call VMCCIS(vec_len,fA,cBha1,mode)
      call VMCCISI(vec_len,fA_I,inca,cBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VCCIS(vec_len,fA,cBha2)
      call VCCISI(vec_len,fA_I,inca,cBha2_I,incb)

      do i=1,10
          if(cBha0(i) .ne. cBha1(i)) then
              print *,"Error! Difference between VCCIS and"
              print *," VMCCIS in VML_HA mode detected"
              stop 1
          endif
          if(cBha0_I((i-1)*incb+1) .ne. cBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCCISI and"
              print *," VMCCISI in VML_HA mode detected"
              stop 1
          endif
          if(cBha1(i) .ne. cBha2(i)) then
              print *,"Error! Difference between VCCIS and"
              print *," VMCCIS in VML_HA mode detected"
              stop 1
          endif
          if(cBha1_I((i-1)*incb+1) .ne. cBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCCISI and"
              print *," VMCCISI in VML_HA mode detected"
              stop 1
          endif
          if(cBla1(i) .ne. cBla2(i)) then
              print *,"Error! Difference between VCCIS and"
              print *," VMCCIS in VML_LA mode detected"
              stop 1
          endif
          if(cBla1_I((i-1)*incb+1) .ne. cBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCCISI and"
              print *," VMCCISI in VML_LA mode detected"
              stop 1
          endif
          if(cBep1(i) .ne. cBep2(i)) then
              print *,"Error! Difference between VCCIS and"
              print *," VMCCIS in VML_EP mode detected"
              stop 1
          endif
          if(cBep1_I((i-1)*incb+1) .ne. cBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCCISI and"
              print *," VMCCISI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vcCIS test/example program"
      print *,""
      print *,"   Argument                     vcCIS"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,fA(i),"    ",cBha0(i)
            CurRMS=crelerr(cB(i),cBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vcCISI test/example program"
      print *,"   Argument                     vcCIS"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,fA_I((i-1)*inca+1),"    ",cBha0_I((i-1)*incb+1)
            CurRMS_I=crelerr(cB_I((i-1)*incb+1),cBha0_I((i-1)*incb+1))
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

10    format(E15.7,A5,E15.7,E15.7)
11    format(A,F25.16)

      end
