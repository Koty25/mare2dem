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
!    vcPowx  Example Program Text
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

      cA( 1)=(0.1000,7.0000)
      cA( 2)=(0.8666,6.2333)
      cA( 3)=(1.6333,5.4666)
      cA( 4)=(2.4000,4.6999)
      cA( 5)=(3.1666,3.9333)
      cA( 6)=(3.9333,3.1666)
      cA( 7)=(4.7000,2.3999)
      cA( 8)=(5.4666,1.6333)
      cA( 9)=(6.2333,0.8666)
      cA(10)=(7.0000,0.0999)
      cB( 1)=(1.4660081731727082e+001,-4.2457664528605692e+000)
      cB( 2)=(6.6546281414811581e+000,-1.3285049956433053e+001)
      cB( 3)=(-7.0387482404509374e+000,-1.4767324238678261e+001)
      cB( 4)=(-2.1117849271454592e+001,-3.4930927491228259e+000)
      cB( 5)=(-2.4424384995282274e+001,2.4394003894154913e+001)
      cB( 6)=(4.1042341669135736e+000,6.8100665478272219e+001)
      cB( 7)=(1.0377464969389693e+002,1.1602925928663136e+002)
      cB( 8)=(3.5030956185511837e+002,1.3713892387072394e+002)
      cB( 9)=(8.9261612350130133e+002,7.6109263759848176e+001)
      cB(10)=(2.0121482556255414e+003,-1.3059791392822208e+002)

      do i=1,10
          cA_I((i-1)*inca+1)=cA(i)
          cB_I((i-1)*incb+1)=cB(i)
      end do

      call VCPOWX(vec_len,cA,cA(6),cBha0)
      call VCPOWXI(vec_len,cA_I,inca,cA(6),cBha0_I,incb)

      mode=VML_EP
      call VMCPOWX(vec_len,cA,cA(6),cBep1,mode)
      call VMCPOWXI(vec_len,cA_I,inca,cA(6),cBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VCPOWX(vec_len,cA,cA(6),cBep2)
      call VCPOWXI(vec_len,cA_I,inca,cA(6),cBep2_I,incb)

      mode=VML_LA
      call VMCPOWX(vec_len,cA,cA(6),cBla1,mode)
      call VMCPOWXI(vec_len,cA_I,inca,cA(6),cBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VCPOWX(vec_len,cA,cA(6),cBla2)
      call VCPOWXI(vec_len,cA_I,inca,cA(6),cBla2_I,incb)

      mode=VML_HA
      call VMCPOWX(vec_len,cA,cA(6),cBha1,mode)
      call VMCPOWXI(vec_len,cA_I,inca,cA(6),cBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VCPOWX(vec_len,cA,cA(6),cBha2)
      call VCPOWXI(vec_len,cA_I,inca,cA(6),cBha2_I,incb)

      do i=1,10
          if(cBha0(i) .ne. cBha1(i)) then
              print *,"Error! Difference between VCPOWX and"
              print *," VMCPOWX in VML_HA mode detected"
              stop 1
          endif
          if(cBha0_I((i-1)*incb+1) .ne. cBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCPOWXI and"
              print *," VMCPOWXI in VML_HA mode detected"
              stop 1
          endif
          if(cBha1(i) .ne. cBha2(i)) then
              print *,"Error! Difference between VCPOWX and"
              print *," VMCPOWX in VML_HA mode detected"
              stop 1
          endif
          if(cBha1_I((i-1)*incb+1) .ne. cBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCPOWXI and"
              print *," VMCPOWXI in VML_HA mode detected"
              stop 1
          endif
          if(cBla1(i) .ne. cBla2(i)) then
              print *,"Error! Difference between VCPOWX and"
              print *," VMCPOWX in VML_LA mode detected"
              stop 1
          endif
          if(cBla1_I((i-1)*incb+1) .ne. cBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCPOWXI and"
              print *," VMCPOWXI in VML_LA mode detected"
              stop 1
          endif
          if(cBep1(i) .ne. cBep2(i)) then
              print *,"Error! Difference between VCPOWX and"
              print *," VMCPOWX in VML_EP mode detected"
              stop 1
          endif
          if(cBep1_I((i-1)*incb+1) .ne. cBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCPOWXI and"
              print *," VMCPOWXI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vcPowx test/example program"
      print *,""
      print *,"           Argument                     vcPowx"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,cA(i),"    ",cBha0(i)
            CurRMS=crelerr(cB(i),cBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vcPowxI test/example program"
      print *,"           Argument                     vcPowx"
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
