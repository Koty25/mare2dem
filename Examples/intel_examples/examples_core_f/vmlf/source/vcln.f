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
!    vcLn  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      complex(kind=4) cA(10)
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
      complex(kind=4) cA_I(10*inca)
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

      cA( 1)=(0.1000,10000.0000)
      cA( 2)=(1111.2000,8888.9004)
      cA( 3)=(2222.2998,7777.7998)
      cA( 4)=(3333.3999,6666.7002)
      cA( 5)=(4444.5000,5555.6001)
      cA( 6)=(5555.6001,4444.5000)
      cA( 7)=(6666.7002,3333.3999)
      cA( 8)=(7777.7998,2222.3000)
      cA( 9)=(8888.9004,1111.2000)
      cA(10)=(10000.0000,0.1000)
      cB( 1)=(9.2103403720261827e+000,1.5707863569259644e+000)
      cB( 2)=(9.1003119337631233e+000,1.4464316368103027e+000)
      cB( 3)=(8.9982670046073512e+000,1.2924882173538208e+000)
      cB( 4)=(8.9164550570889674e+000,1.1071426868438721e+000)
      cB( 5)=(8.8699115947627423e+000,8.9605319499969482e-001)
      cB( 6)=(8.8699115947627423e+000,6.7474311590194702e-001)
      cB( 7)=(8.9164550570889674e+000,4.6365359425544739e-001)
      cB( 8)=(8.9982670128991238e+000,2.7830815315246582e-001)
      cB( 9)=(9.1003119337631233e+000,1.2436468899250031e-001)
      cB(10)=(9.2103403720261827e+000,9.9999997473787516e-006)

      do i=1,10
          cA_I((i-1)*inca+1)=cA(i)
          cB_I((i-1)*incb+1)=cB(i)
      end do

      call VCLN(vec_len,cA,cBha0)
      call VCLNI(vec_len,cA_I,inca,cBha0_I,incb)

      mode=VML_EP
      call VMCLN(vec_len,cA,cBep1,mode)
      call VMCLNI(vec_len,cA_I,inca,cBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VCLN(vec_len,cA,cBep2)
      call VCLNI(vec_len,cA_I,inca,cBep2_I,incb)

      mode=VML_LA
      call VMCLN(vec_len,cA,cBla1,mode)
      call VMCLNI(vec_len,cA_I,inca,cBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VCLN(vec_len,cA,cBla2)
      call VCLNI(vec_len,cA_I,inca,cBla2_I,incb)

      mode=VML_HA
      call VMCLN(vec_len,cA,cBha1,mode)
      call VMCLNI(vec_len,cA_I,inca,cBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VCLN(vec_len,cA,cBha2)
      call VCLNI(vec_len,cA_I,inca,cBha2_I,incb)

      do i=1,10
          if(cBha0(i) .ne. cBha1(i)) then
              print *,"Error! Difference between VCLN and"
              print *," VMCLN in VML_HA mode detected"
              stop 1
          endif
          if(cBha0_I((i-1)*incb+1) .ne. cBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCLNI and"
              print *," VMCLNI in VML_HA mode detected"
              stop 1
          endif
          if(cBha1(i) .ne. cBha2(i)) then
              print *,"Error! Difference between VCLN and"
              print *," VMCLN in VML_HA mode detected"
              stop 1
          endif
          if(cBha1_I((i-1)*incb+1) .ne. cBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCLNI and"
              print *," VMCLNI in VML_HA mode detected"
              stop 1
          endif
          if(cBla1(i) .ne. cBla2(i)) then
              print *,"Error! Difference between VCLN and"
              print *," VMCLN in VML_LA mode detected"
              stop 1
          endif
          if(cBla1_I((i-1)*incb+1) .ne. cBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCLNI and"
              print *," VMCLNI in VML_LA mode detected"
              stop 1
          endif
          if(cBep1(i) .ne. cBep2(i)) then
              print *,"Error! Difference between VCLN and"
              print *," VMCLN in VML_EP mode detected"
              stop 1
          endif
          if(cBep1_I((i-1)*incb+1) .ne. cBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCLNI and"
              print *," VMCLNI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vcLn test/example program"
      print *,""
      print *,"           Argument                     vcLn"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,cA(i),"    ",cBha0(i)
            CurRMS=crelerr(cB(i),cBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vcLnI test/example program"
      print *,"           Argument                     vcLn"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,cA_I((i-1)*inca+1),"    ",cBha0_I((i-1)*incb+1)
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

10    format(E15.7,E15.7,A5,E15.7,E15.7)
11    format(A,F25.16)

      end
