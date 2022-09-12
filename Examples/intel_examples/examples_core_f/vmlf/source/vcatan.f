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
!    vcAtan  Example Program Text
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
      cB( 1)=(-1.5707463267948132e+000,4.9999998736893758e-005)
      cB( 2)=(-1.5707320410808823e+000,6.4285712142009288e-005)
      cB( 3)=(-1.5707063267961685e+000,9.0000001364387572e-005)
      cB( 4)=(-1.5706463267889845e+000,1.4999999257270247e-004)
      cB( 5)=(-1.5703463267231603e+000,4.4999993406236172e-004)
      cB( 6)=(1.5703463267231603e+000,-4.4999993406236172e-004)
      cB( 7)=(1.5706463267889845e+000,-1.4999999257270247e-004)
      cB( 8)=(1.5707063267961685e+000,-9.0000001364387572e-005)
      cB( 9)=(1.5707320410808823e+000,-6.4285712142009288e-005)
      cB(10)=(1.5707463267948132e+000,-4.9999998736893758e-005)

      do i=1,10
          cA_I((i-1)*inca+1)=cA(i)
          cB_I((i-1)*incb+1)=cB(i)
      end do

      call VCATAN(vec_len,cA,cBha0)
      call VCATANI(vec_len,cA_I,inca,cBha0_I,incb)

      mode=VML_EP
      call VMCATAN(vec_len,cA,cBep1,mode)
      call VMCATANI(vec_len,cA_I,inca,cBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VCATAN(vec_len,cA,cBep2)
      call VCATANI(vec_len,cA_I,inca,cBep2_I,incb)

      mode=VML_LA
      call VMCATAN(vec_len,cA,cBla1,mode)
      call VMCATANI(vec_len,cA_I,inca,cBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VCATAN(vec_len,cA,cBla2)
      call VCATANI(vec_len,cA_I,inca,cBla2_I,incb)

      mode=VML_HA
      call VMCATAN(vec_len,cA,cBha1,mode)
      call VMCATANI(vec_len,cA_I,inca,cBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VCATAN(vec_len,cA,cBha2)
      call VCATANI(vec_len,cA_I,inca,cBha2_I,incb)

      do i=1,10
          if(cBha0(i) .ne. cBha1(i)) then
              print *,"Error! Difference between VCATAN and"
              print *," VMCATAN in VML_HA mode detected"
              stop 1
          endif
          if(cBha0_I((i-1)*incb+1) .ne. cBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCATANI and"
              print *," VMCATANI in VML_HA mode detected"
              stop 1
          endif
          if(cBha1(i) .ne. cBha2(i)) then
              print *,"Error! Difference between VCATAN and"
              print *," VMCATAN in VML_HA mode detected"
              stop 1
          endif
          if(cBha1_I((i-1)*incb+1) .ne. cBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCATANI and"
              print *," VMCATANI in VML_HA mode detected"
              stop 1
          endif
          if(cBla1(i) .ne. cBla2(i)) then
              print *,"Error! Difference between VCATAN and"
              print *," VMCATAN in VML_LA mode detected"
              stop 1
          endif
          if(cBla1_I((i-1)*incb+1) .ne. cBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCATANI and"
              print *," VMCATANI in VML_LA mode detected"
              stop 1
          endif
          if(cBep1(i) .ne. cBep2(i)) then
              print *,"Error! Difference between VCATAN and"
              print *," VMCATAN in VML_EP mode detected"
              stop 1
          endif
          if(cBep1_I((i-1)*incb+1) .ne. cBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VCATANI and"
              print *," VMCATANI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vcAtan test/example program"
      print *,""
      print *,"           Argument                     vcAtan"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,cA(i),"    ",cBha0(i)
            CurRMS=crelerr(cB(i),cBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vcAtanI test/example program"
      print *,"           Argument                     vcAtan"
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
