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
!    vcExp  Example Program Text
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

      cA( 1)=(-17.0000,18.0000)
      cA( 2)=(-13.1111,14.1111)
      cA( 3)=(-9.2222,10.2222)
      cA( 4)=(-5.3333,6.3333)
      cA( 5)=(-1.4444,2.4444)
      cA( 6)=(2.4444,-1.4444)
      cA( 7)=(6.3333,-5.3333)
      cA( 8)=(10.2222,-9.2222)
      cA( 9)=(14.1111,-13.1111)
      cA(10)=(18.0000,-17.0000)
      cB( 1)=(2.7336700468037275e-008,-3.1090404917222259e-008)
      cB( 2)=(5.2718026213238074e-008,2.0219670204824070e-006)
      cB( 3)=(-6.9031750320313866e-005,-7.0712150773033500e-005)
      cB( 4)=(4.8220487634171522e-003,2.4185802612919360e-004)
      cB( 5)=(-1.8084268297689621e-001,1.5145584940910339e-001)
      cB( 6)=(1.4526703648953678e+000,-1.1431704521179199e+001)
      cB( 7)=(3.2754686308712672e+002,4.5792468261718750e+002)
      cB( 8)=(-2.6944639799017234e+004,-5.5343007812500000e+003)
      cB( 9)=(1.1494097911982329e+006,-6.9640050000000000e+005)
      cB(10)=(-1.8067216284192696e+007,6.3125328000000000e+007)

      do i=1,10
          cA_I((i-1)*inca+1)=cA(i)
          cB_I((i-1)*incb+1)=cB(i)
      end do

      call VCEXP(vec_len,cA,cBha0)
      call VCEXPI(vec_len,cA_I,inca,cBha0_I,incb)

      mode=VML_EP
      call VMCEXP(vec_len,cA,cBep1,mode)
      call VMCEXPI(vec_len,cA_I,inca,cBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VCEXP(vec_len,cA,cBep2)
      call VCEXPI(vec_len,cA_I,inca,cBep2_I,incb)

      mode=VML_LA
      call VMCEXP(vec_len,cA,cBla1,mode)
      call VMCEXPI(vec_len,cA_I,inca,cBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VCEXP(vec_len,cA,cBla2)
      call VCEXPI(vec_len,cA_I,inca,cBla2_I,incb)

      mode=VML_HA
      call VMCEXP(vec_len,cA,cBha1,mode)
      call VMCEXPI(vec_len,cA_I,inca,cBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VCEXP(vec_len,cA,cBha2)
      call VCEXPI(vec_len,cA_I,inca,cBha2_I,incb)

      do i=1,10
          if(cBha0(i) .ne. cBha1(i)) then
              print *,"Error! Difference between VCEXP and"
              print *," VMCEXP in VML_HA mode detected"
              stop 1
          endif
          if(cBha0_I((i-1)*incb+1) .ne. cBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCEXPI and"
              print *," VMCEXPI in VML_HA mode detected"
              stop 1
          endif
          if(cBha1(i) .ne. cBha2(i)) then
              print *,"Error! Difference between VCEXP and"
              print *," VMCEXP in VML_HA mode detected"
              stop 1
          endif
          if(cBha1_I((i-1)*incb+1) .ne. cBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCEXPI and"
              print *," VMCEXPI in VML_HA mode detected"
              stop 1
          endif
          if(cBla1(i) .ne. cBla2(i)) then
              print *,"Error! Difference between VCEXP and"
              print *," VMCEXP in VML_LA mode detected"
              stop 1
          endif
          if(cBla1_I((i-1)*incb+1) .ne. cBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCEXPI and"
              print *," VMCEXPI in VML_LA mode detected"
              stop 1
          endif
          if(cBep1(i) .ne. cBep2(i)) then
              print *,"Error! Difference between VCEXP and"
              print *," VMCEXP in VML_EP mode detected"
              stop 1
          endif
          if(cBep1_I((i-1)*incb+1) .ne. cBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCEXPI and"
              print *," VMCEXPI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vcExp test/example program"
      print *,""
      print *,"           Argument                     vcExp"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,cA(i),"    ",cBha0(i)
            CurRMS=crelerr(cB(i),cBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vcExpI test/example program"
      print *,"           Argument                     vcExp"
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
