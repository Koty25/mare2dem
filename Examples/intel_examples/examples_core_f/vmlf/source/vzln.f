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
!    vzLn  Example Program Text
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

      zA( 1)=(0.1000d+000,10000.0000d+000)
      zA( 2)=(1111.1999d+000,8888.9000d+000)
      zA( 3)=(2222.2999d+000,7777.8000d+000)
      zA( 4)=(3333.3999d+000,6666.7000d+000)
      zA( 5)=(4444.5000d+000,5555.6000d+000)
      zA( 6)=(5555.6000d+000,4444.5000d+000)
      zA( 7)=(6666.7000d+000,3333.4000d+000)
      zA( 8)=(7777.8000d+000,2222.3000d+000)
      zA( 9)=(8888.9000d+000,1111.2000d+000)
      zA(10)=(10000.0000d+000,0.1000d+000)
      zB( 1)=(9.2103403720261827d+000,1.5707863267948969d+000)
      zB( 2)=(9.1003118897854645d+000,1.4464316510411834d+000)
      zB( 3)=(8.9982670310606725d+000,1.2924881891553042d+000)
      zB( 4)=(8.9164550335110988d+000,1.1071427298419583d+000)
      zB( 5)=(8.8699115840444964d+000,8.9605318947080836d-001)
      zB( 6)=(8.8699115840444964d+000,6.7474313732408819d-001)
      zB( 7)=(8.9164550395111224d+000,4.6365360895280644d-001)
      zB( 8)=(8.9982670344569815d+000,2.7830814952629329d-001)
      zB( 9)=(9.1003118911701844d+000,1.2436468683059552d-001)
      zB(10)=(9.2103403720261827d+000,9.9999999996666679d-006)

      do i=1,10
          zA_I((i-1)*inca+1)=zA(i)
          zB_I((i-1)*incb+1)=zB(i)
      end do

      call VZLN(vec_len,zA,zBha0)
      call VZLNI(vec_len,zA_I,inca,zBha0_I,incb)

      mode=VML_EP
      call VMZLN(vec_len,zA,zBep1,mode)
      call VMZLNI(vec_len,zA_I,inca,zBep1_I,incb,mode)

      tmode=VML_EP
      tmode=VMLSETMODE(tmode)
      call VZLN(vec_len,zA,zBep2)
      call VZLNI(vec_len,zA_I,inca,zBep2_I,incb)

      mode=VML_LA
      call VMZLN(vec_len,zA,zBla1,mode)
      call VMZLNI(vec_len,zA_I,inca,zBla1_I,incb,mode)

      tmode=VML_LA
      tmode=VMLSETMODE(tmode)
      call VZLN(vec_len,zA,zBla2)
      call VZLNI(vec_len,zA_I,inca,zBla2_I,incb)

      mode=VML_HA
      call VMZLN(vec_len,zA,zBha1,mode)
      call VMZLNI(vec_len,zA_I,inca,zBha1_I,incb,mode)

      tmode=VML_HA
      tmode=VMLSETMODE(tmode)
      call VZLN(vec_len,zA,zBha2)
      call VZLNI(vec_len,zA_I,inca,zBha2_I,incb)

      do i=1,10
          if(zBha0(i) .ne. zBha1(i)) then
              print *,"Error! Difference between VZLN and"
              print *," VMZLN in VML_HA mode detected"
              stop 1
          endif
          if(zBha0_I((i-1)*incb+1) .ne. zBha1_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZLNI and"
              print *," VMZLNI in VML_HA mode detected"
              stop 1
          endif
          if(zBha1(i) .ne. zBha2(i)) then
              print *,"Error! Difference between VZLN and"
              print *," VMZLN in VML_HA mode detected"
              stop 1
          endif
          if(zBha1_I((i-1)*incb+1) .ne. zBha2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZLNI and"
              print *," VMZLNI in VML_HA mode detected"
              stop 1
          endif
          if(zBla1(i) .ne. zBla2(i)) then
              print *,"Error! Difference between VZLN and"
              print *," VMZLN in VML_LA mode detected"
              stop 1
          endif
          if(zBla1_I((i-1)*incb+1) .ne. zBla2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZLNI and"
              print *," VMZLNI in VML_LA mode detected"
              stop 1
          endif
          if(zBep1(i) .ne. zBep2(i)) then
              print *,"Error! Difference between VZLN and"
              print *," VMZLN in VML_EP mode detected"
              stop 1
          endif
          if(zBep1_I((i-1)*incb+1) .ne. zBep2_I((i-1)*incb+1)) then
              print *,"Error! Difference between VZLNI and"
              print *," VMZLNI in VML_EP mode detected"
              stop 1
          endif
      end do

      print *,"vzLn test/example program"
      print *,""
      print *,"           Argument                     vzLn"
      print *,"======================================================", &
     &        "========================"
      do i=1,vec_len
            print 10,zA(i),"    ",zBha0(i)
            CurRMS=zrelerr(zB(i),zBha0(i))
            if(CurRMS>MaxRMS) MaxRMS=CurRMS
      end do

      print *,""
      print *,"vzLnI test/example program"
      print *,"           Argument                     vzLn"
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
