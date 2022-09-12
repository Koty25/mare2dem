!===============================================================================
! Copyright 2010-2020 Intel Corporation.
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

      MODULE DF_EXAMPLE_RESCHECK

        REAL(4),PARAMETER ::  EPSILON_SINGLE = 1e-02
        REAL(8),PARAMETER ::  EPSILON_DOUBLE = 1e-05

        INTERFACE
          REAL(4) FUNCTION sComputePPSplineDer(der,n,ny,piece_idx,       &
     &      x,x_i,data_set,sorder,sppcoeff)
            INTEGER,INTENT(IN) :: der
            INTEGER,INTENT(IN) :: n
            INTEGER,INTENT(IN) :: ny
            INTEGER,INTENT(IN) :: piece_idx
            REAL(4),INTENT(IN) :: x
            REAL(4),INTENT(IN) :: x_i
            INTEGER,INTENT(IN) :: data_set
            INTEGER,INTENT(IN) :: sorder
            REAL(4),INTENT(IN) :: sppcoeff((n-1)*ny*sorder)
          END FUNCTION
        END INTERFACE

        INTERFACE
          REAL(8) FUNCTION dComputePPSplineDer(der,n,ny,piece_idx,       &
     &      x,x_i,data_set,sorder,sppcoeff)
            INTEGER,INTENT(IN) :: der
            INTEGER,INTENT(IN) :: n
            INTEGER,INTENT(IN) :: ny
            INTEGER,INTENT(IN) :: piece_idx
            REAL(8),INTENT(IN) :: x
            REAL(8),INTENT(IN) :: x_i
            INTEGER,INTENT(IN) :: data_set
            INTEGER,INTENT(IN) :: sorder
            REAL(8),INTENT(IN) :: sppcoeff((n-1)*ny*sorder)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sFindCells(n,x,nsite,site,bin)
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(4),INTENT(IN) :: site(*)
            INTEGER,INTENT(OUT) :: bin(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dFindCells(n,x,nsite,site,bin)
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(8),INTENT(IN) :: site(*)
            INTEGER,INTENT(OUT) :: bin(*)
          END FUNCTION
        END INTERFACE


        INTERFACE
          INTEGER FUNCTION sCheckQuadBreakPoints(n,x,ny,y,sppcoeff,left, &
     &      right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: y(*)
            REAL(4),INTENT(IN) :: sppcoeff(*)
            REAL(4),INTENT(OUT) :: left(*)
            REAL(4),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckQuadBreakPoints(n,x,ny,y,sppcoeff,left, &
     &      right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: y(*)
            REAL(8),INTENT(IN) :: sppcoeff(*)
            REAL(8),INTENT(OUT) :: left(*)
            REAL(8),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE
        INTERFACE

          INTEGER FUNCTION sCheckSubbQNodes(n,x,nic,ic,ny,sppcoeff,left, &
     &      right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: nic
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: ic(*)
            REAL(4),INTENT(IN) :: sppcoeff(*)
            REAL(4),INTENT(OUT) :: left(*)
            REAL(4),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckSubbQNodes(n,x,nic,ic,ny,sppcoeff,left, &
     &      right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: nic
            REAL(8),INTENT(IN) :: ic(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: sppcoeff(*)
            REAL(8),INTENT(OUT) :: left(*)
            REAL(8),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckSubbBreakPoints(n,x,ny,y,sppcoeff,      &
     &      spline_val)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: y(*)
            REAL(4),INTENT(IN) :: sppcoeff(*)
            REAL(4),INTENT(OUT) :: spline_val(*)
          END FUNCTION
        END INTERFACE

       INTERFACE
          INTEGER FUNCTION dCheckSubbBreakPoints(n,x,ny,y,sppcoeff,      &
     &      spline_val)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: y(*)
            REAL(8),INTENT(IN) :: sppcoeff(*)
            REAL(8),INTENT(OUT) :: spline_val(*)
          END FUNCTION
        END INTERFACE


        INTERFACE
          INTEGER FUNCTION sCheckQuad1stDerConsistency(n,x,ny,sppcoeff,  &
     &      left,right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: sppcoeff(*)
            REAL(4),INTENT(OUT) :: left(*)
            REAL(4),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckQuad1stDerConsistency(n,x,ny,sppcoeff,  &
     &      left,right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: sppcoeff(*)
            REAL(8),INTENT(OUT) :: left(*)
            REAL(8),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckQuadSubb1stDerConsistency(n,x,qnode,ny, &
     &     sppcoeff, left,right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            REAL(4),INTENT(IN) :: qnode(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: sppcoeff(*)
            REAL(4),INTENT(OUT) :: left(*)
            REAL(4),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckQuadSubb1stDerConsistency(n,x,qnode,ny, &
     &     sppcoeff, left,right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            REAL(8),INTENT(IN) :: qnode(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: sppcoeff(*)
            REAL(8),INTENT(OUT) :: left(*)
            REAL(8),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckQuadBC(n,x,ny,sppcoeff,bc_type,bc)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: sppcoeff(*)
            INTEGER,INTENT(IN) :: bc_type
            REAL(4),INTENT(IN) :: bc(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckQuadBC(n,x,ny,sppcoeff,bc_type,bc)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: sppcoeff(*)
            INTEGER,INTENT(IN) :: bc_type
            REAL(8),INTENT(IN) :: bc(*)
          END FUNCTION
        END INTERFACE



        INTERFACE
          INTEGER FUNCTION sCheckCubBC(n,x,ny,sppcoeff,bc_type,bc)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: sppcoeff(*)
            INTEGER,INTENT(IN) :: bc_type
            REAL(4),INTENT(IN),OPTIONAL :: bc(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckCubBC(n,x,ny,sppcoeff,bc_type,bc)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: sppcoeff(*)
            INTEGER,INTENT(IN) :: bc_type
            REAL(8),INTENT(IN),OPTIONAL :: bc(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckCubIntegrRes(n,x,ny,scoeff,nlim,llim,   &
     &          rlim,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN)  :: n
            REAL(4),INTENT(IN)  :: x(*)
            INTEGER,INTENT(IN)  :: ny
            REAL(4),INTENT(IN)  :: scoeff(*)
            INTEGER,INTENT(IN)  :: nlim
            REAL(4),INTENT(IN)  :: llim(*)
            REAL(4),INTENT(IN)  :: rlim(*)
            REAL(4),INTENT(IN)  :: r(*)
            REAL(4),INTENT(OUT) :: ref_r(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckCubIntegrRes(n,x,ny,scoeff,nlim,llim,   &
     &          rlim,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN)  :: n
            REAL(8),INTENT(IN)  :: x(*)
            INTEGER,INTENT(IN)  :: ny
            REAL(8),INTENT(IN)  :: scoeff(*)
            INTEGER,INTENT(IN)  :: nlim
            REAL(8),INTENT(IN)  :: llim(*)
            REAL(8),INTENT(IN)  :: rlim(*)
            REAL(8),INTENT(IN)  :: r(*)
            REAL(8),INTENT(OUT) :: ref_r(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckCubBreakPoints(n,x,ny,y,sppcoeff,left,  &
     &      right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: y(*)
            REAL(4),INTENT(IN) :: sppcoeff(*)
            REAL(4),INTENT(OUT) :: left(*)
            REAL(4),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckCubBreakPoints(n,x,ny,y,sppcoeff,left,  &
     &      right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: y(*)
            REAL(8),INTENT(IN) :: sppcoeff(*)
            REAL(8),INTENT(OUT) :: left(*)
            REAL(8),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckCub1stDerConsistency(n,x,ny,sppcoeff,   &
     &      left,right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: sppcoeff(*)
            REAL(4),INTENT(OUT) :: left(*)
            REAL(4),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckCub1stDerConsistency(n,x,ny,sppcoeff,   &
     &      left,right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: sppcoeff(*)
            REAL(8),INTENT(OUT) :: left(*)
            REAL(8),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckCub2ndDerConsistency(n,x,ny,sppcoeff,   &
     &      left,right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: sppcoeff(*)
            REAL(4),INTENT(OUT) :: left(*)
            REAL(4),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckCub2ndDerConsistency(n,x,ny,sppcoeff,   &
     &      left,right)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: sppcoeff(*)
            REAL(8),INTENT(OUT) :: left(*)
            REAL(8),INTENT(OUT) :: right(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckPPSplineInterpRes(n,x,ny,sorder,scoeff, &
     &      nsite,site,ndorder,dorder,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            INTEGER,INTENT(IN) :: sorder
            REAL(4),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(4),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(4),INTENT(IN) :: r(*)
            REAL(4),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckPPSplineInterpRes(n,x,ny,sorder,scoeff, &
     &      nsite,site,ndorder,dorder,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            INTEGER,INTENT(IN) :: sorder
            REAL(8),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(8),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(8),INTENT(IN) :: r(*)
            REAL(8),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckLinInterpRes(n,x,ny,scoeff,nsite,site,  &
     &      ndorder,dorder,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(4),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(4),INTENT(IN) :: r(*)
            REAL(4),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckLinInterpRes(n,x,ny,scoeff,nsite,site,  &
     &      ndorder,dorder,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(8),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(8),INTENT(IN) :: r(*)
            REAL(8),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckCubInterpRes(n,x,ny,scoeff,nsite,site,  &
     &      ndorder,dorder,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(4),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(4),INTENT(IN) :: r(*)
            REAL(4),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckCubInterpRes(n,x,ny,scoeff,nsite,site,  &
     &      ndorder,dorder,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(8),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(8),INTENT(IN) :: r(*)
            REAL(8),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckQuadInterpRes(n,x,ny,scoeff,nsite,site, &
     &      ndorder,dorder,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(4),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(4),INTENT(IN) :: r(*)
            REAL(4),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckQuadInterpRes(n,x,ny,scoeff,nsite,site, &
     &      ndorder,dorder,r,ref_r)
            USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(8),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(8),INTENT(IN) :: r(*)
            REAL(8),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckLinSplineBreakPoints(n,x,ny,y,sppcoeff)
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: y(*)
            REAL(4),INTENT(IN) :: sppcoeff(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckLinSplineBreakPoints(n,x,ny,y,sppcoeff)
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: y(*)
            REAL(8),INTENT(IN) :: sppcoeff(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION sCheckFifthOrderInterpRes(n,x,ny,scoeff,      &
     &      nsite,site,ndorder,dorder,r,ref_r)
          USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(4),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(4),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(4),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(4),INTENT(IN) :: r(*)
            REAL(4),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER FUNCTION dCheckFifthOrderInterpRes(n,x,ny,scoeff,      &
     &      nsite,site,ndorder,dorder,r,ref_r)
          USE MKL_DF_TYPE
            INTEGER,INTENT(IN) :: n
            REAL(8),INTENT(IN) :: x(*)
            INTEGER,INTENT(IN) :: ny
            REAL(8),INTENT(IN) :: scoeff(*)
            INTEGER,INTENT(IN) :: nsite
            REAL(8),INTENT(IN) :: site(*)
            INTEGER,INTENT(IN) :: ndorder
            INTEGER,INTENT(IN) :: dorder(*)
            REAL(8),INTENT(IN) :: r(*)
            REAL(8),INTENT(OUT) :: ref_r(*)
         END FUNCTION
        END INTERFACE

        INTERFACE
          REAL(4) FUNCTION sComputeCubSplineIntegr(n,ny,piece_idx,llim,  &
     &      rlim,x_i,data_set,sorder,sppcoeff)
            INTEGER,INTENT(IN) :: n
            INTEGER,INTENT(IN) :: ny
            INTEGER,INTENT(IN) :: piece_idx
            REAL(4),INTENT(IN) :: llim
            REAL(4),INTENT(IN) :: rlim
            REAL(4),INTENT(IN) :: x_i
            INTEGER,INTENT(IN) :: data_set
            INTEGER,INTENT(IN) :: sorder
            REAL(4),INTENT(IN) :: sppcoeff((n-1)*ny*sorder)
          END FUNCTION
        END INTERFACE

        INTERFACE
          REAL(8) FUNCTION dComputeCubSplineIntegr(n,ny,piece_idx,llim,  &
     &      rlim,x_i, data_set,sorder,sppcoeff)
            INTEGER,INTENT(IN) :: n
            INTEGER,INTENT(IN) :: ny
            INTEGER,INTENT(IN) :: piece_idx
            REAL(8),INTENT(IN) :: llim
            REAL(8),INTENT(IN) :: rlim
            REAL(8),INTENT(IN) :: x_i
            INTEGER,INTENT(IN) :: data_set
            INTEGER,INTENT(IN) :: sorder
            REAL(8),INTENT(IN) :: sppcoeff((n-1)*ny*sorder)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER(KIND=4) FUNCTION sCheckQuad1stDerBC(n,x,ny,sorder,     &
     &   sppcoeff,bc,isLeft)
             INTEGER,INTENT(IN)  :: n
             INTEGER,INTENT(IN)  :: ny
             INTEGER,INTENT(IN)  :: sorder
             INTEGER,INTENT(IN)  :: isLeft
             REAL(4),INTENT(IN) :: bc
             REAL(4),INTENT(IN) :: x(*)
             REAL(4),INTENT(IN) :: sppcoeff(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER(KIND=4) FUNCTION sCheckQuad2ndDerBC(n,x,ny,sorder,     &
     &   sppcoeff,bc,isLeft)
             INTEGER,INTENT(IN)  :: n
             INTEGER,INTENT(IN)  :: ny
             INTEGER,INTENT(IN)  :: sorder
             INTEGER,INTENT(IN)  :: isLeft
             REAL(4),INTENT(IN) :: bc
             REAL(4),INTENT(IN) :: x(*)
             REAL(4),INTENT(IN) :: sppcoeff(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER(KIND=4) FUNCTION dCheckQuad1stDerBC(n,x,ny,sorder,     &
     &   sppcoeff,bc,isLeft)
             INTEGER,INTENT(IN)  :: n
             INTEGER,INTENT(IN)  :: ny
             INTEGER,INTENT(IN)  :: sorder
             INTEGER,INTENT(IN)  :: isLeft
             REAL(8),INTENT(IN) :: bc
             REAL(8),INTENT(IN) :: x(*)
             REAL(8),INTENT(IN) :: sppcoeff(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER(KIND=4) FUNCTION dCheckQuad2ndDerBC(n,x,ny,sorder,     &
     &   sppcoeff,bc,isLeft)
             INTEGER,INTENT(IN)  :: n
             INTEGER,INTENT(IN)  :: ny
             INTEGER,INTENT(IN)  :: sorder
             INTEGER,INTENT(IN)  :: isLeft
             REAL(8),INTENT(IN) :: bc
             REAL(8),INTENT(IN) :: x(*)
             REAL(8),INTENT(IN) :: sppcoeff(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER(KIND=4) FUNCTION sCheckCub1stDerBC(n,x,ny,sorder,      &
     &   sppcoeff,bc,isLeft)
             INTEGER,INTENT(IN)  :: n
             INTEGER,INTENT(IN)  :: ny
             INTEGER,INTENT(IN)  :: sorder
             INTEGER,INTENT(IN)  :: isLeft
             REAL(4),INTENT(IN) :: bc
             REAL(4),INTENT(IN) :: x(*)
             REAL(4),INTENT(IN) :: sppcoeff(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER(KIND=4) FUNCTION dCheckCub1stDerBC(n,x,ny,sorder,      &
     &   sppcoeff,bc,isLeft)
             INTEGER,INTENT(IN)  :: n
             INTEGER,INTENT(IN)  :: ny
             INTEGER,INTENT(IN)  :: sorder
             INTEGER,INTENT(IN)  :: isLeft
             REAL(8),INTENT(IN) :: bc
             REAL(8),INTENT(IN) :: x(*)
             REAL(8),INTENT(IN) :: sppcoeff(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER(KIND=4) FUNCTION sCheckCub2ndDerBC(n,x,ny,sorder,      &
     &   sppcoeff,bc,isLeft)
             INTEGER,INTENT(IN)  :: n
             INTEGER,INTENT(IN)  :: ny
             INTEGER,INTENT(IN)  :: sorder
             INTEGER,INTENT(IN)  :: isLeft
             REAL(4),INTENT(IN) :: bc
             REAL(4),INTENT(IN) :: x(*)
             REAL(4),INTENT(IN) :: sppcoeff(*)
          END FUNCTION
        END INTERFACE

        INTERFACE
          INTEGER(KIND=4) FUNCTION dCheckCub2ndDerBC(n,x,ny,sorder,      &
     &   sppcoeff,bc,isLeft)
             INTEGER,INTENT(IN) :: n
             REAL(8),INTENT(IN) :: x(*)
             INTEGER,INTENT(IN) :: ny
             INTEGER,INTENT(IN) :: sorder
             REAL(8),INTENT(IN) :: sppcoeff(*)
             REAL(8),INTENT(IN) :: bc
             INTEGER,INTENT(IN) :: isLeft
          END FUNCTION
        END INTERFACE

      END MODULE DF_EXAMPLE_RESCHECK


      ! --------------------------------------------------
      FUNCTION sComputePPSplineDer(der,n,ny,piece_idx,x,x_i,data_set,    &
     & sorder,sppcoeff) RESULT (der_val)
      ! --------------------------------------------------

      ! - - - arg types - - -
        REAL(4) :: der_val
        INTEGER :: der,n,ny,piece_idx,data_set,sorder
        REAL(4),INTENT(IN) :: x,x_i
        REAL(4),TARGET :: sppcoeff((n-1)*ny*sorder)
      ! - - - local declarations - - -
        REAL(4),POINTER :: coeff(:)
        REAL(8) :: t,xpow
        REAL(8) :: d_der_val
        INTEGER :: i,k,mult

      ! - - - begin - - -
        coeff => sppcoeff(data_set * sorder * ( n - 1 ) +                &
     &                    sorder * piece_idx + 1:)
        t = DBLE(x) - DBLE(x_i)
        xpow = 1.0d0
        d_der_val = 0.0d0
        do i = der, sorder-1
          mult = 1
          do k = i-der+1, i
            mult = mult * k
          end do
          d_der_val = d_der_val + DBLE(mult) * DBLE(coeff(i+1)) * xpow
          xpow = xpow * t
        end do

        der_val = REAL(d_der_val,4)

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dComputePPSplineDer(der,n,ny,piece_idx,x,x_i,data_set,    &
     & sorder,sppcoeff)  RESULT (der_val)
      ! --------------------------------------------------

      ! - - - arg types - - -
        REAL(8) :: der_val
        INTEGER :: der,n,ny,piece_idx,data_set,sorder
        REAL(8) :: x,x_i
        REAL(8),TARGET :: sppcoeff((n-1)*ny*sorder)
      ! - - - local declarations - - -
        REAL(8),POINTER :: coeff(:)
        REAL(8) :: t,xpow
        INTEGER :: i,k,mult

      ! - - - begin - - -
        coeff => sppcoeff(data_set * sorder * ( n - 1 ) +                &
     &                    sorder * piece_idx + 1:)
        t = x - x_i
        xpow = 1.0d0
        der_val = 0.0d0
        do i = der, sorder-1
          mult = 1
          do k = i-der+1, i
            mult = mult * k
          end do
          der_val = der_val + DBLE(mult) * coeff(i+1) * xpow
          xpow = xpow * t
        end do

        RETURN

      END FUNCTION


      FUNCTION sCheckLinSplineBreakPoints(n,x,ny,y,sppcoeff)             &
     &  RESULT (errcode)

       USE MKL_DF_TYPE
       USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckLinSplineBreakPoints

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(4) :: x(*),y(*)
        REAL(4),TARGET :: sppcoeff( (n-1)*ny*DF_PP_LINEAR )
      ! - - - local declarations - - -
        INTEGER :: i,j
        INTEGER :: sorder
        REAL(4) :: spline_val_right,spline_val_left,dif
        REAL(4),POINTER :: coeff(:),coeff0(:)

      ! - - - begin - - -
        errcode = 0
        sorder = DF_PP_LINEAR

        DO i = 0, ny-1
          coeff => sppcoeff(i * sorder * ( n - 1 )+1:)
          spline_val_left = ( coeff(1) + coeff(2) * ( x(1) - x(1) ) )
          dif = spline_val_left - y(i*n + 1)
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_SINGLE  ) THEN
            errcode = -1
            RETURN
          END IF
          DO j = 1, n-2
            spline_val_right = ( coeff(sorder*(j-1) + 1)                 &
     &       + coeff(sorder*(j-1) + 2) * ( x(j+1) - x(j) ) )
            spline_val_left  = ( coeff(sorder*j + 1)                     &
     &       + coeff(sorder*j + 2) * ( x(j+1) - x(j+1) ) )
            dif = spline_val_right - spline_val_left
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) THEN
              errcode = -1
              RETURN
            END IF
          END DO
          spline_val_right = coeff(sorder*(n-2) + 1)                     &
     &     + coeff(sorder*(n-2) + 2) * ( x(n) - x(n-1) )
          dif = spline_val_right - y(i*n + n)
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_SINGLE  ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckLinSplineBreakPoints(n,x,ny,y,sppcoeff)             &
     &  RESULT (errcode)
      ! --------------------------------------------------

       USE MKL_DF_TYPE
       USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckLinSplineBreakPoints

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(8) :: x(*),y(*)
        REAL(8),TARGET :: sppcoeff( (n-1)*ny*DF_PP_LINEAR )
      ! - - - local declarations - - -
        INTEGER :: i,j
        INTEGER :: sorder
        REAL(8) :: spline_val_right,spline_val_left,dif
        REAL(8),POINTER :: coeff(:),coeff0(:)

      ! - - - begin - - -
        errcode = 0
        sorder = DF_PP_LINEAR
        DO i = 0, ny-1
          coeff => sppcoeff(i * sorder * ( n - 1 )+1:)
          spline_val_left = ( coeff(1) + coeff(2) * ( x(1) - x(1) ) )
          dif = spline_val_left - y(i*n + 1)
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_DOUBLE  ) THEN
            errcode = -1
            RETURN
          END IF
          DO j = 1, n-2
            spline_val_right = ( coeff(sorder*(j-1) + 1)                 &
     &       + coeff(sorder*(j-1) + 2) * ( x(j+1) - x(j) ) )
            spline_val_left  = ( coeff(sorder*j + 1)                     &
     &       + coeff(sorder*j + 2) * ( x(j+1) - x(j+1) ) )
            dif = spline_val_right - spline_val_left
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE ) THEN
              errcode = -1
              RETURN
            END IF
          END DO
          spline_val_right = coeff(sorder*(n-2) + 1)                     &
     &     + coeff(sorder*(n-2) + 2) * ( x(n) - x(n-1) )
          dif = spline_val_right - y(i*n + n)
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_DOUBLE  ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckCubBreakPoints(n,x,ny,y,sppcoeff,left,right)        &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckCubBreakPoints

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(4) :: x(*),y(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,RowFlag,errnums
        INTEGER :: sorder
        REAL(4) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_CUBIC

        DO i = 0, ny-1
          DO j = 1, n-1
            right(i*(n-1) + j) = sComputePPSplineDer( 0, n, ny, j-1,     &
     &          x(j), x(j), i, sorder, sppcoeff )
            dif = right(i*(n-1) + j) - y(i * n + j+0)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1

            left(i*(n-1) + j) = sComputePPSplineDer( 0, n, ny, j-1,      &
     &          x(j+1), x(j), i, sorder, sppcoeff )
            dif = left(i*(n-1) + j) - y(i * n + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
          END DO
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckCubBreakPoints(n,x,ny,y,sppcoeff,left,right)        &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckCubBreakPoints

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,scoeffhint
        REAL(8) :: x(*),y(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,RowFlag,errnums
        INTEGER :: sorder
        REAL(8) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_CUBIC

        DO i = 0, ny-1
          DO j = 1, n-1
            right(i*(n-1) + j) = dComputePPSplineDer( 0, n, ny, j-1,     &
     &          x(j), x(j), i, sorder, sppcoeff )
            dif = right(i*(n-1) + j) - y(i * n + j+0)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE ) errnums = errnums+1

            left(i*(n-1) + j) = dComputePPSplineDer( 0, n, ny, j-1,      &
     &          x(j+1), x(j), i, sorder, sppcoeff )
            dif = left(i*(n-1) + j) - y(i * n + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE ) errnums = errnums+1
          END DO
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckCub1stDerConsistency(n,x,ny,sppcoeff,left,right)    &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckCub1stDerConsistency

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,scoeffhint
        REAL(4) :: x(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(4) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_CUBIC
        der = 1

        DO i = 0, ny-1
          right(i*(n-1) + 1) = sComputePPSplineDer(der, n, ny, 0, x(1),  &
     &     x(1), i, sorder, sppcoeff)

          DO j = 1, n-2
            left(i*(n-1) + j) = sComputePPSplineDer(der, n, ny, j-1,     &
     &       x(j+1), x(j), i, sorder, sppcoeff)
            right(i*(n-1) + j+1) = sComputePPSplineDer(der, n, ny, j,    &
     &       x(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(n-1) + j) - right(i*(n-1) + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE  ) errnums = errnums+1
          END DO
        left(i*(n-1) + n-1) = sComputePPSplineDer(der, n, ny, n-2, x(n), &
     &   x(n-1), i, sorder, sppcoeff)
        END DO

      IF (errnums > 0) errcode = -1

      RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckCub1stDerConsistency(n,x,ny,sppcoeff,left,right)    &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckCub1stDerConsistency

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,scoeffhint
        REAL(8) :: x(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(8) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_CUBIC
        der = 1

        DO i = 0, ny-1
          right(i*(n-1) + 1) = dComputePPSplineDer(der, n, ny, 0, x(1),  &
     &     x(1), i, sorder, sppcoeff)

          DO j = 1, n-2
            left(i*(n-1) + j) = dComputePPSplineDer(der, n, ny, j-1,     &
     &       x(j+1), x(j), i, sorder, sppcoeff)
            right(i*(n-1) + j+1) = dComputePPSplineDer(der, n, ny, j,    &
     &       x(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(n-1) + j) - right(i*(n-1) + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE  ) errnums = errnums+1
          END DO
        left(i*(n-1) + n-1) = dComputePPSplineDer(der, n, ny, n-2, x(n), &
     &   x(n-1), i, sorder, sppcoeff)
        END DO

      IF (errnums > 0) errcode = -1

      RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckCub2ndDerConsistency(n,x,ny,sppcoeff,left,right)    &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckCub2ndDerConsistency

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(4) :: x(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(4) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_CUBIC
        der = 2

        DO i = 0,ny-1
          right(i*(n-1) + 1) = sComputePPSplineDer(der, n, ny, 0, x(1),  &
     &      x(1), i, sorder, sppcoeff)

          DO j = 1, n-2
            left(i*(n-1) + j) = sComputePPSplineDer(der, n, ny, j-1,     &
     &        x(j+1), x(j), i, sorder, sppcoeff)
            right(i*(n-1) + j+1) = sComputePPSplineDer(der, n, ny, j,    &
     &        x(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(n-1) + j) - right(i*(n-1) + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE  ) errnums = errnums+1
          END DO
          left(i*(n-1) + n-1) = sComputePPSplineDer(der, n, ny, n-2,     &
     &      x(n), x(n-1), i, sorder, sppcoeff)
        END DO

      IF (errnums > 0) errcode = -1

      RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckCub2ndDerConsistency(n,x,ny,sppcoeff,left,right)    &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckCub2ndDerConsistency

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(8) :: x(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(8) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_CUBIC
        der = 2

        DO i = 0,ny-1
          right(i*(n-1) + 1) = dComputePPSplineDer(der, n, ny, 0, x(1),  &
     &      x(1), i, sorder, sppcoeff)

          DO j = 1, n-2
            left(i*(n-1) + j) = dComputePPSplineDer(der, n, ny, j-1,     &
     &        x(j+1), x(j), i, sorder, sppcoeff)
            right(i*(n-1) + j+1) = dComputePPSplineDer(der, n, ny, j,    &
     &        x(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(n-1) + j) - right(i*(n-1) + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE  ) errnums = errnums+1
          END DO
          left(i*(n-1) + n-1) = dComputePPSplineDer(der, n, ny, n-2,     &
     &      x(n), x(n-1), i, sorder, sppcoeff)
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckCubBC(n,x,ny,sppcoeff,bc_type,bc) RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckCubBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,bc_type
        REAL(4) :: x(*),sppcoeff(*),bc(*)
      ! - - - local declarations - - -
        INTEGER :: sorder

        sorder = DF_PP_CUBIC
      ! - - - begin - - -
        errcode = 0
        IF ( bc_type .EQ. DF_BC_NOT_A_KNOT ) THEN
          RETURN
        END IF
        IF ( IAND(bc_type,DF_BC_FREE_END) .EQ. DF_BC_FREE_END ) THEN
          errcode = sCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, 0.0,  &
     &         1 )
          IF ( errcode < 0 ) RETURN
          errcode = sCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, 0.0,  &
     &         0 )
        ELSE
          IF ( IAND(bc_type, DF_BC_1ST_LEFT_DER)                         &
     &     .EQ. DF_BC_1ST_LEFT_DER) THEN
            errcode = sCheckCub1stDerBC( n, x, ny, sorder, sppcoeff,     &
     &          bc(1), 1 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_1ST_RIGHT_DER)                        &
     &     .EQ. DF_BC_1ST_RIGHT_DER ) THEN
            errcode = sCheckCub1stDerBC( n, x, ny, sorder, sppcoeff,     &
     &          bc(2), 0 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_2ND_LEFT_DER)                         &
     &     .EQ. DF_BC_2ND_LEFT_DER ) THEN
            errcode = sCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff,     &
     &          bc(1), 1 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_2ND_RIGHT_DER)                        &
     &     .EQ. DF_BC_2ND_RIGHT_DER ) THEN
            errcode = sCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff,     &
     &           bc(2), 0 )
            IF ( errcode < 0 ) RETURN
          END IF
        END IF

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckCubBC(n,x,ny,sppcoeff,bc_type,bc) RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckCubBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,bc_type
        REAL(8) :: x(*),sppcoeff(*),bc(*)
      ! - - - local declarations - - -
        INTEGER :: sorder

        sorder = DF_PP_CUBIC

      ! - - - begin - - -
        errcode = 0
        IF ( bc_type .EQ. DF_BC_NOT_A_KNOT ) THEN
          RETURN
        END IF
        IF ( IAND(bc_type,DF_BC_FREE_END) .EQ. DF_BC_FREE_END ) THEN
          errcode = dCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff,       &
     &         0.0d0, 1 )
          IF ( errcode < 0 ) RETURN
          errcode = dCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff,       &
     &         0.0d0, 0 )
        ELSE
          IF ( IAND(bc_type, DF_BC_1ST_LEFT_DER)                         &
     &     .EQ. DF_BC_1ST_LEFT_DER) THEN
            errcode = dCheckCub1stDerBC( n, x, ny, sorder, sppcoeff,     &
     &          bc(1), 1 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_1ST_RIGHT_DER)                        &
     &     .EQ. DF_BC_1ST_RIGHT_DER ) THEN
            errcode = dCheckCub1stDerBC( n, x, ny, sorder, sppcoeff,     &
     &          bc(2), 0 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_2ND_LEFT_DER)                         &
     &     .EQ. DF_BC_2ND_LEFT_DER ) THEN
            errcode = dCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff,     &
     &          bc(1), 1 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_2ND_RIGHT_DER)                        &
     &     .EQ. DF_BC_2ND_RIGHT_DER ) THEN
            errcode = dCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff,     &
     &           bc(2), 0 )
            IF ( errcode < 0 ) RETURN
          END IF
        END IF

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sComputeCubSplineIntegr(n,ny,piece_idx,llim,rlim,x_i,     &
     & data_set,sorder,sppcoeff)  RESULT (integr_val)
      ! --------------------------------------------------

      ! - - - arg types - - -
        REAL(4) :: integr_val
        INTEGER,INTENT(IN) :: n,ny,piece_idx,data_set,sorder
        REAL(4),INTENT(IN) :: llim,rlim,x_i
        REAL(4),TARGET,INTENT(IN) :: sppcoeff((n-1)*ny*sorder)
        ! - - - local declarations - - -
        REAL(4) :: coeff0,coeff1,coeff2,coeff3,t
        REAL(4),POINTER :: coeff(:)

        ! - - - begin - - -
        coeff => sppcoeff(data_set * sorder * ( n - 1 )+1:)

        coeff0 = coeff(sorder*piece_idx + 1)
        coeff1 = coeff(sorder*piece_idx + 2)
        coeff2 = coeff(sorder*piece_idx + 3)
        coeff3 = coeff(sorder*piece_idx + 4)

        ! Compute integartion coefficients
        coeff1 = coeff1 * 0.5
        coeff2 = coeff2 / 3.0
        coeff3 = coeff3 * 0.25

        ! Integartion result computation
        t = rlim - x_i
        integr_val = t*(coeff0 + t*(coeff1 + t*(coeff2 + t*coeff3)))

        t = llim - x_i
        integr_val = integr_val -                                        &
     &               t*(coeff0 + t*(coeff1 + t*(coeff2 + t*coeff3)))
        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION dComputeCubSplineIntegr(n,ny,piece_idx,llim,rlim,x_i,     &
     & data_set,sorder,sppcoeff)  RESULT (integr_val)
      ! --------------------------------------------------

      ! - - - arg types - - -
        REAL(8) :: integr_val
        INTEGER,INTENT(IN) :: n,ny,piece_idx,data_set,sorder
        REAL(8),INTENT(IN) :: llim,rlim,x_i
        REAL(8),TARGET,INTENT(IN) :: sppcoeff((n-1)*ny*sorder)
        ! - - - local declarations - - -
        REAL(8) :: coeff0,coeff1,coeff2,coeff3,t
        REAL(8),POINTER :: coeff(:)

        ! - - - begin - - -
        coeff => sppcoeff(data_set * sorder * ( n - 1 )+1:)

        coeff0 = coeff(sorder*piece_idx + 1)
        coeff1 = coeff(sorder*piece_idx + 2)
        coeff2 = coeff(sorder*piece_idx + 3)
        coeff3 = coeff(sorder*piece_idx + 4)

        ! Compute integartion coefficients
        coeff1 = coeff1 * 0.5
        coeff2 = coeff2 / 3.0
        coeff3 = coeff3 * 0.25

        ! Integartion result computation
        t = rlim - x_i
        integr_val = t*(coeff0 + t*(coeff1 + t*(coeff2 + t*coeff3)))

        t = llim - x_i
        integr_val = integr_val -                                        &
     &               t*(coeff0 + t*(coeff1 + t*(coeff2 + t*coeff3)))
        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckCub1stDerBC(n,x,ny,sorder,sppcoeff,bc,isLeft)       &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckCub1stDerBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,isLeft
        REAL(4) :: bc,x(*)
        REAL(4), target :: sppcoeff((n-1)*ny*sorder)

      ! - - - local declarations - - -
        INTEGER :: i,piece_idx,x_idx,der
        REAL(4) :: spline_der_val,dif

      ! - - - begin - - -
        errcode = 0
        der = 1
        IF (isLeft .NE. 0) THEN
          piece_idx = 0
        ELSE
          piece_idx = (n-2)
        END IF
        IF (isLeft .NE. 0) THEN
          x_idx = 0
        ELSE
          x_idx = (n-1)
        END IF
        DO i = 0, ny-1
          spline_der_val = sComputePPSplineDer(der, n, ny, piece_idx,    &
     &     x(x_idx+1), x(piece_idx+1), i, sorder, sppcoeff)
          dif = spline_der_val - bc
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_SINGLE  ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckCub1stDerBC(n,x,ny,sorder,sppcoeff,bc,isLeft)       &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckCub1stDerBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,isLeft
        REAL(8) :: bc,x(*)
        REAL(8), target :: sppcoeff((n-1)*ny*sorder)

      ! - - - local declarations - - -
        INTEGER :: i,piece_idx,x_idx,der
        REAL(8) :: spline_der_val,dif

      ! - - - begin - - -
        errcode = 0
        der = 1
        IF (isLeft .NE. 0) THEN
          piece_idx = 0
        ELSE
          piece_idx = (n-2)
        END IF
        IF (isLeft .NE. 0) THEN
          x_idx = 0
        ELSE
          x_idx = (n-1)
        END IF
        DO i = 0,ny-1
          spline_der_val = dComputePPSplineDer(der, n, ny, piece_idx,    &
     &     x(x_idx+1), x(piece_idx+1), i, sorder, sppcoeff)
          dif = spline_der_val - bc
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_DOUBLE  ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckCub2ndDerBC(n,x,ny,sorder,sppcoeff,bc,isLeft)       &
     &  RESULT(errcode)
      ! --------------------------------------------------
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckCub2ndDerBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,isLeft
        REAL(4) :: bc,x(*)
        REAL(4), target :: sppcoeff((n-1)*ny*sorder)

      ! - - - local declarations - - -
        INTEGER :: i,piece_idx,x_idx,der
        REAL(4) :: spline_der_val,dif

      ! - - - begin - - -
        errcode = 0
        der = 2
        IF (isLeft .NE. 0) THEN
          piece_idx = 0
        ELSE
          piece_idx = (n-2)
        END IF
        IF (isLeft .NE. 0) THEN
          x_idx = 0
        ELSE
          x_idx = (n-1)
        END IF
        DO i = 0,ny-1
          spline_der_val = sComputePPSplineDer(der, n, ny, piece_idx,    &
     &     x(x_idx+1), x(piece_idx+1), i, sorder, sppcoeff)
          dif = spline_der_val - bc
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_SINGLE  ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckCub2ndDerBC(n,x,ny,sorder,sppcoeff,bc,isLeft)       &
     &  RESULT(errcode)
      ! --------------------------------------------------
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckCub2ndDerBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,isLeft
        REAL(8) :: bc,x(*)
        REAL(8), target :: sppcoeff((n-1)*ny*sorder)

      ! - - - local declarations - - -
        INTEGER :: i,piece_idx,x_idx,der
        REAL(8) :: spline_der_val,dif

      ! - - - begin - - -
        errcode = 0
        der = 2
        IF (isLeft .NE. 0) THEN
          piece_idx = 0
        ELSE
          piece_idx = (n-2)
        END IF

        IF (isLeft .NE. 0) THEN
          x_idx = 0
        ELSE
          x_idx = (n-1)
        END IF

        DO i = 0,ny-1
          spline_der_val = dComputePPSplineDer(der, n, ny, piece_idx,    &
     &     x(x_idx+1), x(piece_idx+1), i, sorder, sppcoeff)
          dif = spline_der_val - bc
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_DOUBLE  ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckQuad1stDerBC(n,x,ny,sorder,sppcoeff,bc,isLeft)      &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckQuad1stDerBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,isLeft
        REAL(4) :: bc,x(*),sppcoeff((n-1)*ny*sorder)
      ! - - - local declarations - - -
        INTEGER :: i,piece_idx,x_idx,der
        REAL(4) :: spline_der_val,dif

      ! - - - begin - - -
        errcode = 0
        IF (isLeft .NE. 0) THEN
          piece_idx = 0
          x_idx = 0
        ELSE
          piece_idx = (n-1)
          x_idx = (n-1)
        END IF

        der = 1

        DO i = 0, ny-1
          spline_der_val = sComputePPSplineDer(der, n,ny, piece_idx,     &
     &       x(x_idx+1), x(piece_idx+1), i, sorder, sppcoeff)

          dif = spline_der_val - bc
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_SINGLE  ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION dCheckQuad1stDerBC(n,x,ny,sorder,sppcoeff,bc,isLeft)      &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckQuad1stDerBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,isLeft
        REAL(8) :: bc,x(*),sppcoeff((n-1)*ny*sorder)
      ! - - - local declarations - - -
        INTEGER :: i,piece_idx,x_idx,der
        REAL(8) :: spline_der_val,dif

      ! - - - begin - - -
        errcode = 0
        IF (isLeft .NE. 0) THEN
          piece_idx = 0
          x_idx = 0
        ELSE
          piece_idx = (n-1)
          x_idx = (n-1)
        END IF

        der = 1

        DO i = 0,ny-1
          spline_der_val = dComputePPSplineDer(der, n,ny, piece_idx,     &
     &       x(x_idx+1), x(piece_idx+1), i, sorder, sppcoeff)

          dif = spline_der_val - bc
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_DOUBLE  ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION sCheckQuadBreakPoints(n,x,ny,y,sppcoeff,left,right)       &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckQuadBreakPoints

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(4) :: x(*),y(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(4) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 0

        DO i = 0, ny-1
          DO j = 1, n-1
            right(i*(n-1) + j) = sComputePPSplineDer( der, n, ny, j-1,   &
     &          x(j), x(j), i, sorder, sppcoeff )
            dif = right(i*(n-1) + j) - y(i * n + j)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
            left(i*(n-1) + j) = sComputePPSplineDer( der, n, ny, j-1,    &
     &          x(j+1), x(j), i, sorder, sppcoeff )
            dif = left(i*(n-1) + j) - y(i * n + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
          END DO
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckQuadBreakPoints(n,x,ny,y,sppcoeff,left,right)       &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckQuadBreakPoints

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(8) :: x(*),y(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder, der
        REAL(8) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 0

        DO i = 0, ny-1
          DO j = 1, n-1
            right(i*(n-1) + j) = dComputePPSplineDer( der, n, ny, j-1,   &
     &          x(j), x(j), i, sorder, sppcoeff )
            dif = right(i*(n-1) + j) - y(i * n + j)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE ) errnums = errnums+1
            left(i*(n-1) + j) = dComputePPSplineDer( der, n, ny, j-1,    &
     &          x(j+1), x(j), i, sorder, sppcoeff )
            dif = left(i*(n-1) + j) - y(i * n + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE ) errnums = errnums+1
          END DO
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckSubbQNodes(nx,x,nqnode,qnode,ny,sppcoeff,left,      &
     &  right)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckSubbQNodes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: nx,nqnode,ny
        REAL(4) :: x(*),qnode(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(4) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 0

        DO i = 0,ny-1
          right(i*(nqnode-1) + 1) = sComputePPSplineDer(der, nqnode, ny, &
     &         0, qnode(1), x(1), i, sorder, sppcoeff)

          DO j = 1, nqnode-2
            left(i*(nqnode-1) + j) = sComputePPSplineDer(der, nqnode,    &
     &          ny, j-1, qnode(j+1), x(j), i, sorder, sppcoeff)
            right(i*(nqnode-1) + j+1) = sComputePPSplineDer(der, nqnode, &
     &          ny, j, qnode(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(nqnode-1) + j) - right(i*(nqnode-1)+ j+1)
            IF ( dif < 0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
          END DO
          left(i*(nqnode-1) + nqnode-1) = sComputePPSplineDer(der,       &
     &       nqnode, ny, nqnode-2, qnode(nqnode), x(nqnode-1), i,        &
     &       sorder, sppcoeff)
        END DO

      IF (errnums > 0) errcode = -1

      RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION dCheckSubbQNodes(nx,x,nqnode,qnode,ny,sppcoeff,left,      &
     &  right)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckSubbQNodes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: nx,nqnode,ny
        REAL(8) :: x(*),qnode(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(8) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 0

        DO i = 0,ny-1
          right(i*(nqnode-1) + 1) = dComputePPSplineDer(der, nqnode, ny, &
     &         0, qnode(1), x(1), i, sorder, sppcoeff)

          DO j = 1, nqnode-2
            left(i*(nqnode-1) + j) = dComputePPSplineDer(der, nqnode,    &
     &          ny, j-1, qnode(j+1), x(j), i, sorder, sppcoeff)
            right(i*(nqnode-1) + j+1) = dComputePPSplineDer(der, nqnode, &
     &          ny, j, qnode(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(nqnode-1) + j) - right(i*(nqnode-1)+ j+1)
            IF ( dif < 0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE ) errnums = errnums+1
          END DO
          left(i*(nqnode-1) + nqnode-1) = dComputePPSplineDer(der,       &
     &       nqnode, ny, nqnode-2, qnode(nqnode), x(nqnode-1), i,        &
     &       sorder, sppcoeff)
        END DO

      IF (errnums > 0) errcode = -1

      RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckSubbBreakPoints(n,x,ny,y,sppcoeff,spline_val)       &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckSubbBreakPoints

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(4) :: x(*),y(*),sppcoeff(*),spline_val(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums,sorder,der
        REAL(4) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 0

        DO i = 0,ny-1
          DO j = 0,n-1
            spline_val(i*n + j+1) = sComputePPSplineDer( der, n+1, ny,   &
     &           j, x(j+1), x(j+1), i, sorder, sppcoeff )
            dif = spline_val(i*n + j+1) - y(i*n + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
          END DO
        END DO
        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckSubbBreakPoints(n,x,ny,y,sppcoeff,spline_val)       &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckSubbBreakPoints

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(8) :: x(*),y(*),sppcoeff(*),spline_val(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums,sorder,der
        REAL(8) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 0

        DO i = 0,ny-1
          DO j = 0,n-1
            spline_val(i*n + j+1) = dComputePPSplineDer( der, n+1, ny,   &
     &           j, x(j+1), x(j+1), i, sorder, sppcoeff )
            dif = spline_val(i*n + j+1) - y(i*n + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE ) errnums = errnums+1
          END DO
        END DO
        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckQuad1stDerConsistency(n,x,ny,sppcoeff,left,right)   &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK,ignore_me=>sCheckQuad1stDerConsistency

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(4) :: x(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(4) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 1

        DO i = 0, ny-1
          right(i*(n-1) + 1) = sComputePPSplineDer(der, n, ny, 0, x(1),  &
     &        x(1), i, sorder, sppcoeff)
          DO j = 1, n-2
            left(i*(n-1) + j) = sComputePPSplineDer(der, n, ny, j-1,     &
     &         x(j+1), x(j), i, sorder, sppcoeff)
            right(i*(n-1) + j+1) = sComputePPSplineDer(der, n, ny, j,    &
     &         x(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(n-1) + j) - right(i*(n-1) + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
          END DO
          left(i*(n-1) + n-1) = sComputePPSplineDer(der, n, ny, n-2,     &
     &        x(n), x(n-1), i, sorder, sppcoeff)
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckQuad1stDerConsistency(n,x,ny,sppcoeff,left,right)   &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK,ignore_me=>dCheckQuad1stDerConsistency

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(8) :: x(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(8) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 1

        DO i = 0, ny-1
          right(i*(n-1) + 1) = dComputePPSplineDer(der, n, ny, 0, x(1),  &
     &        x(1), i, sorder, sppcoeff)
          DO j = 1, n-2
            left(i*(n-1) + j) = dComputePPSplineDer(der, n, ny, j-1,     &
     &         x(j+1), x(j), i, sorder, sppcoeff)
            right(i*(n-1) + j+1) = dComputePPSplineDer(der, n, ny, j,    &
     &         x(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(n-1) + j) - right(i*(n-1) + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
          END DO
          left(i*(n-1) + n-1) = dComputePPSplineDer(der, n, ny, n-2,     &
     &        x(n), x(n-1), i, sorder, sppcoeff)
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckQuadSubb1stDerConsistency(n,x,qnode,ny,sppcoeff,    &
     & left,right)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK,                                         &
     & ignore_me=>sCheckQuadSubb1stDerConsistency

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(4) :: x(*),qnode(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(4) :: dif

      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 1

        DO i = 0, ny-1
          right(i*(n-1) + 1) = sComputePPSplineDer(der, n, ny, 0,        &
     &        qnode(1) ,x(1), i, sorder, sppcoeff)
          DO j = 1, n-2
            left (i*(n-1) + j)   = sComputePPSplineDer(der, n, ny, j-1,  &
     &          qnode(j+1), x(j), i, sorder, sppcoeff)
            right(i*(n-1) + j+1) = sComputePPSplineDer(der, n, ny, j,    &
     &         qnode(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(n-1) + j) - right(i*(n-1) + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
          END DO
          left(i*(n-1) + n-1) = sComputePPSplineDer(der, n, ny, n-2,     &
     &        qnode(n), x(n-1), i, sorder, sppcoeff)
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION dCheckQuadSubb1stDerConsistency(n,x,qnode,ny,sppcoeff,    &
     & left,right)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK,                                         &
     & ignore_me => dCheckQuadSubb1stDerConsistency

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny
        REAL(8) :: x(*),qnode(*),sppcoeff(*),left(*),right(*)
      ! - - - local declarations - - -
        INTEGER :: i,j,errnums
        INTEGER :: sorder,der
        REAL(8) :: dif


      ! - - - begin - - -
        errcode = 0
        errnums = 0
        sorder = DF_PP_QUADRATIC
        der = 1

        DO i = 0, ny-1
          right(i*(n-1) + 1) = dComputePPSplineDer(der, n, ny, 0,        &
     &        qnode(1) ,x(1), i, sorder, sppcoeff)
          DO j = 1, n-2
            left (i*(n-1) + j)   = dComputePPSplineDer(der, n, ny, j-1,  &
     &          qnode(j+1), x(j), i, sorder, sppcoeff)
            right(i*(n-1) + j+1) = dComputePPSplineDer(der, n, ny, j,    &
     &         qnode(j+1), x(j+1), i, sorder, sppcoeff)
            dif = left(i*(n-1) + j) - right(i*(n-1) + j+1)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
          END DO
          left(i*(n-1) + n-1) = dComputePPSplineDer(der, n, ny, n-2,     &
     &        qnode(n), x(n-1), i, sorder, sppcoeff)
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION sCheckQuadBC(n,x,ny,sppcoeff,bc_type,bc)                  &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckQuadBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,bc_type
        REAL(4) :: x(*),sppcoeff(*),bc(*)
      ! - - - local declarations - - -
        INTEGER :: i,sorder,der
        REAL(4) :: left_val,spline_left_val,dif

        sorder = DF_PP_QUADRATIC

      ! - - - begin - - -
        errcode = 0
        der = 0

        IF ( IAND(bc_type, DF_BC_Q_VAL) .EQ. DF_BC_Q_VAL ) THEN
          left_val = ( x(1) + x(2) ) / 2.0
          DO i = 0, ny-1
            spline_left_val = sComputePPSplineDer(der, n, ny, 0,         &
     &            left_val, x(1), i, sorder, sppcoeff)
            dif = spline_left_val - bc(1)
            IF (dif < 0.0) dif = -dif
            IF (dif > EPSILON_SINGLE ) errcode = -1
          END DO
        ELSE
          IF ( IAND(bc_type, DF_BC_1ST_LEFT_DER)                         &
     &     .EQ. DF_BC_1ST_LEFT_DER ) THEN
            errcode = sCheckQuad1stDerBC( n, x, ny, sorder, sppcoeff,    &
     &          bc(1), 1 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_1ST_RIGHT_DER)                        &
     &     .EQ. DF_BC_1ST_RIGHT_DER ) THEN

            errcode = sCheckQuad1stDerBC( n, x, ny, sorder, sppcoeff,    &
     &           bc(2), 0 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_2ND_LEFT_DER)                         &
     &     .EQ. DF_BC_2ND_LEFT_DER ) THEN
            errcode = sCheckQuad2ndDerBC( n, x, ny, sorder, sppcoeff,    &
     &           bc(1), 1 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_2ND_RIGHT_DER)                        &
     &     .EQ. DF_BC_2ND_RIGHT_DER ) THEN
            errcode = sCheckQuad2ndDerBC( n, x, ny, sorder, sppcoeff,    &
     &           bc(2), 0 )
            IF ( errcode < 0 ) RETURN
          END IF
        END IF

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION dCheckQuadBC(n,x,ny,sppcoeff,bc_type,bc)                  &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckQuadBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,bc_type
        REAL(8) :: x(*),sppcoeff(*),bc(*)
      ! - - - local declarations - - -
        INTEGER :: i,sorder,der
        REAL(8) :: left_val,spline_left_val,dif

        sorder = DF_PP_QUADRATIC

      ! - - - begin - - -
        errcode = 0
        der = 0

        IF ( IAND(bc_type, DF_BC_Q_VAL) .EQ. DF_BC_Q_VAL ) THEN
          left_val = ( x(1) + x(2) ) / 2.0d0
          DO i = 0, ny-1
            spline_left_val = dComputePPSplineDer(der, n, ny, 0,         &
     &            left_val, x(1), i, sorder, sppcoeff)
            dif = spline_left_val - bc(1)
            IF (dif .LT. 0.0) dif = -dif
            IF (dif .GT. EPSILON_DOUBLE) errcode = -1
          END DO
        ELSE
          IF ( IAND(bc_type, DF_BC_1ST_LEFT_DER)                         &
     &     .EQ. DF_BC_1ST_LEFT_DER ) THEN
            errcode = dCheckQuad1stDerBC( n, x, ny, sorder, sppcoeff,    &
     &           bc(1), 1 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_1ST_RIGHT_DER)                        &
     &     .EQ. DF_BC_1ST_RIGHT_DER ) THEN
            errcode = dCheckQuad1stDerBC( n, x, ny, sorder, sppcoeff,    &
     &           bc(2), 0 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_2ND_LEFT_DER)                         &
     &     .EQ. DF_BC_2ND_LEFT_DER ) THEN
            errcode = dCheckQuad2ndDerBC( n, x, ny, sorder, sppcoeff,    &
     &           bc(1), 1 )
            IF ( errcode < 0 ) RETURN
          END IF
          IF ( IAND(bc_type, DF_BC_2ND_RIGHT_DER)                        &
     &     .EQ. DF_BC_2ND_RIGHT_DER ) THEN
            errcode = dCheckQuad2ndDerBC( n, x, ny, sorder, sppcoeff,    &
     &           bc(2), 0 )
            IF ( errcode < 0 ) RETURN
          END IF
        END IF

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION sCheckQuad2ndDerBC(n,x,ny,sorder,sppcoeff,bc,isLeft)      &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckQuad2ndDerBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,isLeft
        REAL(4) :: bc,x(*),sppcoeff(*)
      ! - - - local declarations - - -
        INTEGER :: i,piece_idx,x_idx,der
        REAL(4) :: spline_der_val,dif

      ! - - - begin - - -
        errcode = 0
        IF (isLeft .NE. 0) THEN
          piece_idx = 0
        ELSE
          piece_idx = (n-1)
        END IF
        IF (isLeft .NE. 0) THEN
          x_idx = 0
        ELSE
          x_idx = (n-1)
        END IF

        der = 2

        DO i = 0,ny-1
          spline_der_val = sComputePPSplineDer(der, n, ny, piece_idx,    &
     &         x(x_idx+1), x(x_idx+1), i, sorder, sppcoeff)
          dif = spline_der_val - bc
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_SINGLE ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckQuad2ndDerBC(n,x,ny,sorder,sppcoeff,bc,isLeft)      &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckQuad2ndDerBC

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,isLeft
        REAL(8) :: bc,x(*),sppcoeff(*)
      ! - - - local declarations - - -
        INTEGER :: i,piece_idx,x_idx,der
        REAL(8) :: spline_der_val,dif

      ! - - - begin - - -
        errcode = 0
        IF (isLeft .NE. 0) THEN
          piece_idx = 0
        ELSE
          piece_idx = (n-1)
        END IF
        IF (isLeft .NE. 0) THEN
          x_idx = 0
        ELSE
          x_idx = (n-1)
        END IF

        der = 2

        DO i = 0,ny-1
          spline_der_val = dComputePPSplineDer(der, n, ny, piece_idx,    &
     &          x(x_idx+1), x(x_idx+1), i, sorder, sppcoeff)
          dif = spline_der_val - bc
          IF ( dif < 0.0 ) dif = -dif
          IF ( dif > EPSILON_DOUBLE ) THEN
            errcode = -1
            RETURN
          END IF
        END DO

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sFindCells(n,x,nsite,site,bin)  RESULT (errcode)
      ! --------------------------------------------------

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,nsite,bin(*)
        REAL(4) :: x(*),site(*)
      ! - - - local declarations - - -
        INTEGER :: i,j

      ! - - - begin - - -
        DO j = 1, nsite
          IF ( site(j) .EQ. x(n) ) THEN
            bin(j) = n
          ELSE
          DO i = 1, n
              bin(j) = n+1
            IF ( site(j) < x(i) ) THEN
              bin(j) = i
              EXIT
            END IF
          END DO
          END IF
        END DO

        errcode = 0
        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dFindCells(n,x,nsite,site,bin)  RESULT (errcode)
      ! --------------------------------------------------

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,nsite,bin(*)
        REAL(8) :: x(*),site(*)
      ! - - - local declarations - - -
        INTEGER :: i, j

      ! - - - begin - - -
        DO j = 1, nsite
          IF ( site(j) .EQ. x(n) ) THEN
            bin(j) = n
          ELSE
          DO i = 1, n
              bin(j) = n+1
              IF ( site(j) < x(i) ) THEN
              bin(j) = i
              EXIT
            END IF
          END DO
          END IF
        END DO

        errcode = 0
        RETURN
      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckPPSplineInterpRes(n,x,ny,sorder,scoeff,nsite,site,  &
     & ndorder,dorder,r,ref_r)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckPPSplineInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,nsite,ndorder,dorder(*)
        REAL(4) :: x(*),scoeff(*),site(*)
        REAL(4),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: i,j,j1,jder,k,errnums,nderorder
        INTEGER :: bin(nsite)
        REAL(4) :: dif
        REAL(4),POINTER :: ref_r_ptr(:),r_ptr(:)

      ! - - - begin - - -
        errcode = 0
        errnums = 0

        errcode = sFindCells(n, x, nsite, site, bin)
        DO i = 1, nsite
          bin(i) = bin(i) - 1
          IF (bin(i) .GT. n-1) bin(i) = n-1
        END DO

        nderorder = 0
        DO i = 1, ndorder
          IF ( dorder(i) .NE. 0 ) nderorder = nderorder+1
        END DO

        DO i = 1, nderorder*nsite*ny
          ref_r(i) = 0
        END DO

        DO i = 0, ny-1
          ref_r_ptr => ref_r(i*nderorder*nsite+1:)
          r_ptr => r(i*nderorder*nsite+1:)
          DO j = 1, nsite
            j1 = 1
            IF ((site(j) .GE. x(1)) .AND. (site(j) .LE. x(n))) THEN
              DO jder = 0, nderorder
                IF ((ndorder > jder) .AND. (dorder(jder+1) .NE. 0)) THEN
                  ref_r_ptr((j-1)*nderorder+j1)=                         &
     &             sComputePPSplineDer(jder,n,                           &
     &             ny,bin(j)-1,site(j), x(bin(j)), i, sorder, scoeff)
                  dif = (ref_r_ptr((j-1)*nderorder+j1) -                 &
     &                   r_ptr    ((j-1)*nderorder+j1))
                  IF ( dif < 0.0 ) dif = -dif
                  IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
                  j1 = j1+1
                END IF
              END DO
            END IF
          END DO
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckPPSplineInterpRes(n,x,ny,sorder,scoeff,nsite,site,  &
     & ndorder,dorder,r,ref_r)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckPPSplineInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,sorder,nsite,ndorder,dorder(*)
        REAL(8) :: x(*),scoeff(*),site(*)
        REAL(8),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: i,j,j1,jder,k,errnums,nderorder
        INTEGER :: bin(nsite)
        REAL(8) :: dif
        REAL(8),POINTER :: ref_r_ptr(:),r_ptr(:)

      ! - - - begin - - -
        errcode = 0
        errnums = 0

        errcode = dFindCells(n, x, nsite, site, bin)
        DO i = 1, nsite
          bin(i) = bin(i) - 1
          IF (bin(i) .GT. n-1) bin(i) = n-1
        END DO

        nderorder = 0
        DO i = 1, ndorder
          IF ( dorder(i) .NE. 0 ) nderorder = nderorder+1
        END DO

        DO i = 1, nderorder*nsite*ny
          ref_r(i) = 0
        END DO

        DO i = 0, ny-1
          ref_r_ptr => ref_r(i*nderorder*nsite+1:)
          r_ptr => r(i*nderorder*nsite+1:)
          DO j = 1, nsite
            j1 = 1
            IF ((site(j) .GE. x(1)) .AND. (site(j) .LE. x(n))) THEN
              DO jder = 0, nderorder
                IF ((ndorder > jder) .AND. (dorder(jder+1) .NE. 0)) THEN
                  ref_r_ptr((j-1)*nderorder+j1)=                         &
     &             dComputePPSplineDer(jder,n,                           &
     &             ny,bin(j)-1,site(j), x(bin(j)), i, sorder, scoeff)
                  dif = (ref_r_ptr((j-1)*nderorder+j1) -                 &
     &                   r_ptr    ((j-1)*nderorder+j1))
                  IF ( dif < 0.0 ) dif = -dif
                  IF ( dif > EPSILON_SINGLE ) errnums = errnums+1
                  j1 = j1+1
                END IF
              END DO
            END IF
          END DO
        END DO

        IF (errnums > 0) errcode = -1

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION sCheckLinInterpRes(n,x,ny,scoeff,nsite,site,              &
     & ndorder,dorder,r,ref_r)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckLinInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,nsite,ndorder,dorder(*)
        REAL(4) :: x(*),scoeff(*),site(*)
        REAL(4),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: sorder

      ! - - - begin - - -
        sorder = DF_PP_LINEAR
        errcode = sCheckPPSplineInterpRes(n, x, ny, sorder, scoeff,      &
     &      nsite, site, ndorder, dorder, r, ref_r)

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckLinInterpRes(n,x,ny,scoeff,nsite,site,              &
     & ndorder,dorder,r,ref_r)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckLinInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,nsite,ndorder,dorder(ndorder)
        REAL(8) :: x(n),scoeff(*),site(nsite)
        REAL(8),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: sorder

      ! - - - begin - - -
        sorder = DF_PP_LINEAR
        errcode = dCheckPPSplineInterpRes(n, x, ny, sorder, scoeff,      &
     &      nsite, site, ndorder, dorder, r, ref_r)

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION sCheckCubInterpRes(n,x,ny,scoeff,nsite,site,              &
     & ndorder,dorder,r,ref_r)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckCubInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,nsite,ndorder,dorder(*)
        REAL(4) :: x(*),scoeff(*),site(*)
        REAL(4),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: sorder

      ! - - - begin - - -
        sorder = DF_PP_CUBIC
        errcode = sCheckPPSplineInterpRes(n, x, ny, sorder, scoeff,      &
     &      nsite, site, ndorder, dorder, r, ref_r)

        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckCubInterpRes(n,x,ny,scoeff,nsite,site,              &
     & ndorder,dorder,r,ref_r)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckCubInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,nsite,ndorder,dorder(ndorder)
        REAL(8) :: x(n),scoeff(*),site(nsite)
        REAL(8),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: sorder

      ! - - - begin - - -
        sorder = DF_PP_CUBIC
        errcode = dCheckPPSplineInterpRes(n, x, ny, sorder, scoeff,      &
     &      nsite, site, ndorder, dorder, r, ref_r)

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION sCheckQuadInterpRes(n,x,ny,scoeff,nsite,site,             &
     & ndorder,dorder,r,ref_r) RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckQuadInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,nsite,ndorder,dorder(ndorder)
        REAL(4) :: x(n),scoeff(*),site(nsite)
        REAL(4),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: sorder

      ! - - - begin - - -
        sorder = DF_PP_QUADRATIC
        errcode = sCheckPPSplineInterpRes(n, x, ny, sorder, scoeff,      &
     &      nsite, site, ndorder, dorder, r, ref_r)

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION dCheckQuadInterpRes(n,x,ny,scoeff,nsite,site,             &
     & ndorder,dorder,r,ref_r) RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckQuadInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,nsite,ndorder,dorder(ndorder)
        REAL(8) :: x(n),scoeff(*),site(nsite)
        REAL(8),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: sorder

      ! - - - begin - - -
        sorder = DF_PP_QUADRATIC
        errcode = dCheckPPSplineInterpRes(n, x, ny, sorder, scoeff,      &
     &      nsite, site, ndorder, dorder, r, ref_r)

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION sCheckCubIntegrRes(n,x,ny,scoeff,nlim,llim,rlim,r,ref_r)  &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckCubIntegrRes

        ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n, ny, nlim
        REAL(4) :: x(n), scoeff(*)
        REAL(4) :: llim(nlim), rlim(nlim)
        REAL(4),TARGET :: r(ny*nlim)
        REAL(4),TARGET :: ref_r(ny*nlim)
      ! - - - local declarations - - -
        INTEGER :: i, j, jbin, tmp_bin, errnums
        INTEGER :: sorder
        INTEGER :: lbin(nlim), rbin(nlim)
        INTEGER :: lbin1, rbin1
        REAL(4) :: sign_flag, tmp_lim, ref_r_tmp, dif
        REAL(4),POINTER :: ref_r_ptr(:),r_ptr(:)

        errcode = 0
        errnums = 0
        sorder = DF_PP_CUBIC

        errcode = sFindCells( n, x, nlim, llim, lbin )
        errcode = sFindCells( n, x, nlim, rlim, rbin )

        DO j = 1, nlim
          lbin(j) = lbin(j) - 1
          rbin(j) = rbin(j) - 1
        END DO

        DO i = 1, ny*nlim
          ref_r(i) = 0.0d0
        END DO

        DO i = 0, ny-1
          ref_r_ptr => ref_r(i*nderorder*nsite+1:)
          r_ptr => r(i*nderorder*nsite+1:)

          DO j = 1, nlim
            sign_flag = 1.0d0
            IF ( rlim(j) < llim(j) ) THEN
              ! Here if right integration limit is lower than left
              tmp_lim = llim(j)
              llim(j) = rlim(j)
              rlim(j) = tmp_lim;

              tmp_bin = lbin(j)
              lbin(j) = rbin(j)
              rbin(j) = tmp_bin

              sign_flag = -1.0;
            END IF

            IF ( lbin(j) .EQ. rbin(j) ) THEN
              ! Here if there are no break points between integration
              ! limits

              ref_r_ptr(j) = sComputeCubSplineIntegr( n, ny, lbin(j)-1,  &
     &              llim(j), rlim(j), x(lbin(j)), i, sorder, scoeff )
            ELSE
              ! Here if there is one or more break points between integration
              ! limits

              ref_r_ptr(j) = sComputeCubSplineIntegr( n, ny, lbin(j)-1,  &
     &              llim(j), x(lbin(j)+1), x(lbin(j)), i, sorder,        &
     &              scoeff )
              DO jbin = lbin(j) + 1, rbin(j)-1

                ref_r_tmp = sComputeCubSplineIntegr( n, ny, jbin-1,      &
     &                x(jbin), x(jbin+1), x(jbin), i, sorder, scoeff )
                ref_r_ptr(j) = ref_r_ptr(j) + ref_r_tmp
              END DO

              ref_r_tmp = sComputeCubSplineIntegr( n, ny, rbin(j)-1,     &
     &              x(rbin(j)), rlim(j), x(rbin(j)), i, sorder,          &
     &              scoeff )
              ref_r_ptr(j) = ref_r_ptr(j) + ref_r_tmp
            END IF

            ref_r_ptr(j) = sign_flag * ref_r_ptr(j)

            dif = ref_r_ptr(j) - r_ptr(j)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_SINGLE ) errnums = errnums + 1
          END DO
        END DO

        IF ( errnums > 0 ) errcode = -1

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION dCheckCubIntegrRes(n,x,ny,scoeff,nlim,llim,rlim,r,ref_r)  &
     &  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckCubIntegrRes

        ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n, ny, nlim
        REAL(8) :: x(n), scoeff(*)
        REAL(8) :: llim(nlim), rlim(nlim)
        REAL(8),TARGET :: r(ny*nlim)
        REAL(8),TARGET :: ref_r(ny*nlim)
      ! - - - local declarations - - -
        INTEGER :: i, j, jbin, tmp_bin, errnums
        INTEGER :: sorder
        INTEGER :: lbin(nlim), rbin(nlim)
        INTEGER :: lbin1, rbin1
        REAL(8) :: sign_flag, tmp_lim, ref_r_tmp, dif
        REAL(8),POINTER :: ref_r_ptr(:),r_ptr(:)

        errcode = 0
        errnums = 0
        sorder = DF_PP_CUBIC

        errcode = dFindCells( n, x, nlim, llim, lbin )
        errcode = dFindCells( n, x, nlim, rlim, rbin )

        DO j = 1, nlim
          lbin(j) = lbin(j) - 1
          rbin(j) = rbin(j) - 1
        END DO

        DO i = 1, ny*nlim
          ref_r(i) = 0.0d0
        END DO

        DO i = 0, ny-1
          ref_r_ptr => ref_r(i*nderorder*nsite+1:)
          r_ptr => r(i*nderorder*nsite+1:)

          DO j = 1, nlim
            sign_flag = 1.0d0
            IF ( rlim(j) < llim(j) ) THEN
              ! Here if right integration limit is lower than left
              tmp_lim = llim(j)
              llim(j) = rlim(j)
              rlim(j) = tmp_lim;

              tmp_bin = lbin(j)
              lbin(j) = rbin(j)
              rbin(j) = tmp_bin

              sign_flag = -1.0;
            END IF

            IF ( lbin(j) .EQ. rbin(j) ) THEN
              ! Here if there are no break points between integration
              ! limits

              ref_r_ptr(j) = dComputeCubSplineIntegr( n, ny, lbin(j)-1,  &
     &              llim(j), rlim(j), x(lbin(j)), i, sorder, scoeff )
            ELSE
              ! Here if there is one or more break points between integration
              ! limits

              ref_r_ptr(j) = dComputeCubSplineIntegr( n, ny, lbin(j)-1,  &
     &              llim(j), x(lbin(j)+1), x(lbin(j)), i, sorder,        &
     &              scoeff )
              DO jbin = lbin(j) + 1, rbin(j)-1

                ref_r_tmp = dComputeCubSplineIntegr( n, ny, jbin-1,      &
     &                x(jbin), x(jbin+1), x(jbin), i, sorder, scoeff )
                ref_r_ptr(j) = ref_r_ptr(j) + ref_r_tmp
              END DO

              ref_r_tmp = dComputeCubSplineIntegr( n, ny, rbin(j)-1,     &
     &              x(rbin(j)), rlim(j), x(rbin(j)), i, sorder,          &
     &              scoeff )
              ref_r_ptr(j) = ref_r_ptr(j) + ref_r_tmp
            END IF

            ref_r_ptr(j) = sign_flag * ref_r_ptr(j)

            dif = ref_r_ptr(j) - r_ptr(j)
            IF ( dif < 0.0 ) dif = -dif
            IF ( dif > EPSILON_DOUBLE ) errnums = errnums + 1
          END DO
        END DO

        IF ( errnums > 0 ) errcode = -1

        RETURN

      END FUNCTION


      ! --------------------------------------------------
      FUNCTION sCheckFifthOrderInterpRes(n,x,ny,scoeff,nsite,site,       &
     & ndorder,dorder,r,ref_r)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => sCheckFifthOrderInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,nsite,ndorder,dorder(*)
        REAL(4) :: x(*),scoeff(*),site(*)
        REAL(4),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: sorder

      ! - - - begin - - -
        sorder = 5
        errcode = sCheckPPSplineInterpRes(n, x, ny, sorder, scoeff,      &
     &      nsite, site, ndorder, dorder, r, ref_r)


        RETURN

      END FUNCTION

      ! --------------------------------------------------
      FUNCTION dCheckFifthOrderInterpRes(n,x,ny,scoeff,nsite,site,       &
     & ndorder,dorder,r,ref_r)  RESULT (errcode)
      ! --------------------------------------------------
        USE MKL_DF_TYPE
        USE DF_EXAMPLE_RESCHECK, ignore_me => dCheckFifthOrderInterpRes

      ! - - - arg types - - -
        INTEGER :: errcode
        INTEGER :: n,ny,nsite,ndorder,dorder(ndorder)
        REAL(8) :: x(n),scoeff(*),site(nsite)
        REAL(8),target :: r(ndorder*nsite*ny),ref_r(ndorder*nsite*ny)
      ! - - - local declarations - - -
        INTEGER :: sorder

      ! - - - begin - - -
        sorder = 5
        errcode = dCheckPPSplineInterpRes(n, x, ny, sorder, scoeff,      &
     &      nsite, site, ndorder, dorder, r, ref_r)

        RETURN

      END FUNCTION
