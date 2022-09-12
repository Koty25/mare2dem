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
!    vcLog10  Example Program Text
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
      cB( 1)=(4.0000000000217151e+000,6.8218386173248291e-001)
      cB( 2)=(3.9522152564316348e+000,6.2817728519439697e-001)
      cB( 3)=(3.9078977067930758e+000,5.6132048368453979e-001)
      cB( 4)=(3.8723672294320828e+000,4.8082599043846130e-001)
      cB( 5)=(3.8521536605751314e+000,3.8915094733238220e-001)
      cB( 6)=(3.8521536605751314e+000,2.9303723573684692e-001)
      cB( 7)=(3.8723672294320828e+000,2.0136220753192902e-001)
      cB( 8)=(3.9078977103941459e+000,1.2086769193410873e-001)
      cB( 9)=(3.9522152564316348e+000,5.4010897874832153e-002)
      cB(10)=(4.0000000000217151e+000,4.3429449760878924e-006)

      do i=1,10
          cA_I((i-1)*inca+1)=cA(i)
          cB_I((i-1)*incb+1)=cB(i)
      end do

      call VCLOG10(vec_len,cA,cBha0)
      call VCLOG10I(vec_len,cA_I,inca,cBha0_I,incb)

      mode=VML_EP
      call VMCLOG10(vec_len,cA,cBep1,mode)
      call VMCLOG10I(vec_len,cA_I,inca,cBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VCLOG10(vec_len,cA,cBep2)
      call VCLOG10I(vec_len,cA_I,inca,cBep2_I,incb)

      mode=VML_LA
      call VMCLOG10(vec_len,cA,cBla1,mode)
      call VMCLOG10I(vec_len,cA_I,inca,cBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VCLOG10(vec_len,cA,cBla2)
      call VCLOG10I(vec_len,cA_I,inca,cBla2_I,incb)

      mode=VML_HA
      call VMCLOG10(vec_len,cA,cBha1,mode)
      call VMCLOG10I(vec_len,cA_I,inca,cBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VCLOG10(vec_len,cA,cBha2)
      call VCLOG10I(vec_len,cA_I,inca,cBha2_I,incb)

      do i=1,10
          if(cBha0(i) .ne. cBha1(i)) then
              print *,"Error! Difference between VCLOG10 and"
              print *," VMCLOG10 in VML_HA mode detected"
              stop 1
          endif
          if(cBha0_I((i-1)*incb+1) .ne. cBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCLOG10I and"
              print *," VMCLOG10I in VML_HA mode detected"
              stop 1
          endif
          if(cBha1(i) .ne. cBha2(i)) then
              print *,"Error! Difference between VCLOG10 and"
              print *," VMCLOG10 in VML_HA mode detected"
              stop 1
          endif
          if(cBha1_I((i-1)*incb+1) .ne. cBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCLOG10I and"
              print *," VMCLOG10I in VML_HA mode detected"
              stop 1
          endif
          if(cBla1(i) .ne. cBla2(i)) then
              print *,"Error! Difference between VCLOG10 and"
              print *," VMCLOG10 in VML_LA mode detected"
              stop 1
          endif
          if(cBla1_I((i-1)*incb+1) .ne. cBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCLOG10I and"
              print *," VMCLOG10I in VML_LA mode detected"
              stop 1
          endif
          if(cBep1(i) .ne. cBep2(i)) then
              print *,"Error! Difference between VCLOG10 and"
              print *," VMCLOG10 in VML_EP mode detected"
              stop 1
          endif
          if(cBep1_I((i-1)*incb+1) .ne. cBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCLOG10I and"
              print *," VMCLOG10I in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vcLog10 test/example program"
      print *,""
      print *,"           Argument                     vcLog10"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,cA(i),"    ",cBha0(i)
            CurRMS=crelerr(cB(i),cBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vcLog10I test/example program"
      print *,"           Argument                     vcLog10"
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
