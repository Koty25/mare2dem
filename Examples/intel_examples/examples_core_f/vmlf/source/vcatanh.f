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
!    vcAtanh  Example Program Text
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

      cA( 1)=(0.7500,0.9000)
      cA( 2)=(0.7666,0.8833)
      cA( 3)=(0.7833,0.8666)
      cA( 4)=(0.8000,0.8500)
      cA( 5)=(0.8166,0.8333)
      cA( 6)=(0.8333,0.8166)
      cA( 7)=(0.8500,0.8000)
      cA( 8)=(0.8666,0.7833)
      cA( 9)=(0.8833,0.7666)
      cA(10)=(0.9000,0.7500)
      cB( 1)=(3.7257323907975143e-001,8.8743013143539429e-001)
      cB( 2)=(3.8548667482649812e-001,8.8805592060089111e-001)
      cB( 3)=(3.9865922844057755e-001,8.8905519247055054e-001)
      cB( 4)=(4.1200697096525207e-001,8.9044255018234253e-001)
      cB( 5)=(4.2558004817630984e-001,8.9211910963058472e-001)
      cB( 6)=(4.3940256872758254e-001,8.9423400163650513e-001)
      cB( 7)=(4.5338331929830056e-001,8.9679872989654541e-001)
      cB( 8)=(4.6758174531596247e-001,8.9971578121185303e-001)
      cB( 9)=(4.8200889281307913e-001,9.0314799547195435e-001)
      cB(10)=(4.9656447251629632e-001,9.0710288286209106e-001)

      do i=1,10
          cA_I((i-1)*inca+1)=cA(i)
          cB_I((i-1)*incb+1)=cB(i)
      end do

      call VCATANH(vec_len,cA,cBha0)
      call VCATANHI(vec_len,cA_I,inca,cBha0_I,incb)

      mode=VML_EP
      call VMCATANH(vec_len,cA,cBep1,mode)
      call VMCATANHI(vec_len,cA_I,inca,cBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VCATANH(vec_len,cA,cBep2)
      call VCATANHI(vec_len,cA_I,inca,cBep2_I,incb)

      mode=VML_LA
      call VMCATANH(vec_len,cA,cBla1,mode)
      call VMCATANHI(vec_len,cA_I,inca,cBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VCATANH(vec_len,cA,cBla2)
      call VCATANHI(vec_len,cA_I,inca,cBla2_I,incb)

      mode=VML_HA
      call VMCATANH(vec_len,cA,cBha1,mode)
      call VMCATANHI(vec_len,cA_I,inca,cBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VCATANH(vec_len,cA,cBha2)
      call VCATANHI(vec_len,cA_I,inca,cBha2_I,incb)

      do i=1,10
          if(cBha0(i) .ne. cBha1(i)) then
              print *,"Error! Difference between VCATANH and"
              print *," VMCATANH in VML_HA mode detected"
              stop 1
          endif
          if(cBha0_I((i-1)*incb+1) .ne. cBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCATANHI and"
              print *," VMCATANHI in VML_HA mode detected"
              stop 1
          endif
          if(cBha1(i) .ne. cBha2(i)) then
              print *,"Error! Difference between VCATANH and"
              print *," VMCATANH in VML_HA mode detected"
              stop 1
          endif
          if(cBha1_I((i-1)*incb+1) .ne. cBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCATANHI and"
              print *," VMCATANHI in VML_HA mode detected"
              stop 1
          endif
          if(cBla1(i) .ne. cBla2(i)) then
              print *,"Error! Difference between VCATANH and"
              print *," VMCATANH in VML_LA mode detected"
              stop 1
          endif
          if(cBla1_I((i-1)*incb+1) .ne. cBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCATANHI and"
              print *," VMCATANHI in VML_LA mode detected"
              stop 1
          endif
          if(cBep1(i) .ne. cBep2(i)) then
              print *,"Error! Difference between VCATANH and"
              print *," VMCATANH in VML_EP mode detected"
              stop 1
          endif
          if(cBep1_I((i-1)*incb+1) .ne. cBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCATANHI and"
              print *," VMCATANHI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vcAtanh test/example program"
      print *,""
      print *,"           Argument                     vcAtanh"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,cA(i),"    ",cBha0(i)
            CurRMS=crelerr(cB(i),cBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vcAtanhI test/example program"
      print *,"           Argument                     vcAtanh"
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
