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
!    vzSqrt  Example Program Text
!*******************************************************************************

      include "_rms.fi"

      program MKL_VML_TEST

      include "mkl_vml.f90"
      real(kind=4) :: srelerr
      real(kind=8) :: drelerr
      real(kind=4) :: crelerr
      real(kind=8) :: zrelerr

      complex(kind=8) zA(10)
      complex(kind=8) zB(10)
      complex(kind=8) zBha0(10)
      complex(kind=8) zBha1(10)
      complex(kind=8) zBha2(10)
      complex(kind=8) zBla1(10)
      complex(kind=8) zBla2(10)
      complex(kind=8) zBep1(10)
      complex(kind=8) zBep2(10)

      real(kind=8) CurRMS,MaxRMS

      integer, parameter :: inca=3
      integer, parameter :: incb=5
      complex(kind=8) zA_I(10*inca)
      complex(kind=8) zB_I(10*incb)
      complex(kind=8) zBha0_I(10*incb)
      complex(kind=8) zBha1_I(10*incb)
      complex(kind=8) zBha2_I(10*incb)
      complex(kind=8) zBla1_I(10*incb)
      complex(kind=8) zBla2_I(10*incb)
      complex(kind=8) zBep1_I(10*incb)
      complex(kind=8) zBep2_I(10*incb)
      real(kind=8) CurRMS_I,MaxRMS_I

      integer(kind=8) mode
      integer tmode
      integer i, vec_len

      vec_len=10
      MaxRMS=0.0
      MaxRMS_I=0.0

      zA( 1)=(0.0000d+000,10000.0000d+000)
      zA( 2)=(1111.1111d+000,8888.8888d+000)
      zA( 3)=(2222.2222d+000,7777.7777d+000)
      zA( 4)=(3333.3333d+000,6666.6666d+000)
      zA( 5)=(4444.4444d+000,5555.5555d+000)
      zA( 6)=(5555.5555d+000,4444.4444d+000)
      zA( 7)=(6666.6666d+000,3333.3333d+000)
      zA( 8)=(7777.7777d+000,2222.2222d+000)
      zA( 9)=(8888.8888d+000,1111.1111d+000)
      zA(10)=(10000.0000d+000,0.0000d+000)
      zB( 1)=(7.0710678118654755d+001,7.0710678118654755d+001)
      zB( 2)=(7.0954827796266002d+001,6.2637660297921109d+001)
      zB( 3)=(7.1802622191669840d+001,5.4160819358644105d+001)
      zB( 4)=(7.3440088338943667d+001,4.5388470730262000d+001)
      zB( 5)=(7.6023111008727710d+001,3.6538596134024317d+001)
      zB( 6)=(7.9593146422574208d+001,2.7919768219763874d+001)
      zB( 7)=(8.4024479916461544d+001,1.9835489034350779d+001)
      zB( 8)=(8.9069603701822302d+001,1.2474638415588544d+001)
      zB( 9)=(9.4464154246982829d+001,5.8811255383440271d+000)
      zB(10)=(1.0000000000000000d+002,0.0000000000000000d+000)

      do i=1,10
          zA_I((i-1)*inca+1)=zA(i)
          zB_I((i-1)*incb+1)=zB(i)
      end do

      call VZSQRT(vec_len,zA,zBha0)
      call VZSQRTI(vec_len,zA_I,inca,zBha0_I,incb)

      mode=VML_EP
      call VMZSQRT(vec_len,zA,zBep1,mode)
      call VMZSQRTI(vec_len,zA_I,inca,zBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZSQRT(vec_len,zA,zBep2)
      call VZSQRTI(vec_len,zA_I,inca,zBep2_I,incb)

      mode=VML_LA
      call VMZSQRT(vec_len,zA,zBla1,mode)
      call VMZSQRTI(vec_len,zA_I,inca,zBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZSQRT(vec_len,zA,zBla2)
      call VZSQRTI(vec_len,zA_I,inca,zBla2_I,incb)

      mode=VML_HA
      call VMZSQRT(vec_len,zA,zBha1,mode)
      call VMZSQRTI(vec_len,zA_I,inca,zBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZSQRT(vec_len,zA,zBha2)
      call VZSQRTI(vec_len,zA_I,inca,zBha2_I,incb)

      do i=1,10
          if(zBha0(i) .ne. zBha1(i)) then
              print *,"Error! Difference between VZSQRT and"
              print *," VMZSQRT in VML_HA mode detected"
              stop 1
          endif
          if(zBha0_I((i-1)*incb+1) .ne. zBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZSQRTI and"
              print *," VMZSQRTI in VML_HA mode detected"
              stop 1
          endif
          if(zBha1(i) .ne. zBha2(i)) then
              print *,"Error! Difference between VZSQRT and"
              print *," VMZSQRT in VML_HA mode detected"
              stop 1
          endif
          if(zBha1_I((i-1)*incb+1) .ne. zBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZSQRTI and"
              print *," VMZSQRTI in VML_HA mode detected"
              stop 1
          endif
          if(zBla1(i) .ne. zBla2(i)) then
              print *,"Error! Difference between VZSQRT and"
              print *," VMZSQRT in VML_LA mode detected"
              stop 1
          endif
          if(zBla1_I((i-1)*incb+1) .ne. zBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZSQRTI and"
              print *," VMZSQRTI in VML_LA mode detected"
              stop 1
          endif
          if(zBep1(i) .ne. zBep2(i)) then
              print *,"Error! Difference between VZSQRT and"
              print *," VMZSQRT in VML_EP mode detected"
              stop 1
          endif
          if(zBep1_I((i-1)*incb+1) .ne. zBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZSQRTI and"
              print *," VMZSQRTI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzSqrt test/example program"
      print *,""
      print *,"           Argument                     vzSqrt"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,zA(i),"    ",zBha0(i)
            CurRMS=zrelerr(zB(i),zBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzSqrtI test/example program"
      print *,"           Argument                     vzSqrt"
      print *,"======================================================", &
     &        "========================"
      do i=1,10
            print 10,zA_I((i-1)*inca+1),"    ",zBha0_I((i-1)*incb+1)
            CurRMS_I=zrelerr(zB_I((i-1)*incb+1),zBha0_I((i-1)*incb+1))
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

10    format(E15.7,E15.7,A5,E15.7,E15.7)
11    format(A,F25.16)

      end
