!===============================================================================
! Copyright 2014-2020 Intel Corporation.
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

!
!     Calculate double precision GaussianMV sample characteristics
!
      SUBROUTINE dCalculateGaussianMVSampleCharacteristics( ndim, n, r,  &
     &      dbS, dbS2, dbMean, dbVar, dbCovXY, dbCovXZ, dbCovYZ )
        integer,intent(in)          :: ndim
        integer,intent(in)          :: n
        real(kind=8),intent(in)     :: r(ndim,n)
        real(kind=8),intent(out)    :: dbS(ndim)
        real(kind=8),intent(out)    :: dbS2(ndim)
        real(kind=8),intent(out)    :: dbMean(ndim)
        real(kind=8),intent(out)    :: dbVar(ndim)
        real(kind=8),intent(out)    :: dbCovXY
        real(kind=8),intent(out)    :: dbCovXZ
        real(kind=8),intent(out)    :: dbCovYZ

        integer i, j
        real(kind=8) dbSXY, dbSXZ, dbSYZ
        real(kind=8) dn

        dn = n
        do i = 1, ndim
          dbS(i) = 0.0D0
          dbS2(i) = 0.0D0
        end do

        dbSXY = 0.0D0
        dbSXZ = 0.0D0
        dbSYZ = 0.0D0

        do i = 1, n
          do j = 1, ndim
            dbS(j)  = dbS(j)  + r(j, i)
            dbS2(j) = dbS2(j) + r(j, i) * r(j, i)
          end do
          dbSXY = dbSXY + r(1, i) * r(2, i)
          dbSXZ = dbSXZ + r(1, i) * r(3, i)
          dbSYZ = dbSYZ + r(2, i) * r(3, i)
        end do

        do j = 1, ndim
          dbMean(j) = dbS(j)/dn
          dbVar(j) = dbS2(j)/dn - dbMean(j)*dbMean(j)
        end do

        dbCovXY = dbSXY/dn - dbMean(1)*dbMean(2)
        dbCovXZ = dbSXZ/dn - dbMean(1)*dbMean(3)
        dbCovYZ = dbSYZ/dn - dbMean(2)*dbMean(3)
      END SUBROUTINE

!
!     Calculate single precision GaussianMV sample characteristics
!
      SUBROUTINE sCalculateGaussianMVSampleCharacteristics( ndim, n, r,  &
     &      dbS, dbS2, dbMean, dbVar, dbCovXY, dbCovXZ, dbCovYZ )
        integer,intent(in)          :: ndim
        integer,intent(in)          :: n
        real(kind=4),intent(in)     :: r(ndim,n)
        real(kind=4),intent(out)    :: dbS(ndim)
        real(kind=4),intent(out)    :: dbS2(ndim)
        real(kind=4),intent(out)    :: dbMean(ndim)
        real(kind=4),intent(out)    :: dbVar(ndim)
        real(kind=4),intent(out)    :: dbCovXY
        real(kind=4),intent(out)    :: dbCovXZ
        real(kind=4),intent(out)    :: dbCovYZ

        integer i, j
        real(kind=4) dbSXY, dbSXZ, dbSYZ
        real(kind=4) dn

        dn = n
        do i = 1, ndim
          dbS(i) = 0.0
          dbS2(i) = 0.0
        end do

        dbSXY = 0.0
        dbSXZ = 0.0
        dbSYZ = 0.0

        do i = 1, n
          do j = 1, ndim
            dbS(j)  = dbS(j)  + r(j, i)
            dbS2(j) = dbS2(j) + r(j, i) * r(j, i)
          end do
          dbSXY = dbSXY + r(1, i) * r(2, i)
          dbSXZ = dbSXZ + r(1, i) * r(3, i)
          dbSYZ = dbSYZ + r(2, i) * r(3, i)
        end do

        do j = 1, ndim
          dbMean(j) = dbS(j)/dn
          dbVar(j) = dbS2(j)/dn - dbMean(j)*dbMean(j)
        end do

        dbCovXY = dbSXY/dn - dbMean(1)*dbMean(2)
        dbCovXZ = dbSXZ/dn - dbMean(1)*dbMean(3)
        dbCovYZ = dbSYZ/dn - dbMean(2)*dbMean(3)
      END SUBROUTINE

!
!     Check that double precision GaussianMV sample characteristics agree with theory
!
      FUNCTION dGaussianMVCheckResults( ndim, n, a, c, dbMean, dbVar, S, &
     &      D2, Q, DeltaM, DeltaD )
        integer,intent(in)          :: ndim
        integer,intent(in)          :: n
        real(kind=8),intent(in)     :: a(ndim)
        real(kind=8),intent(in)     :: c(ndim,ndim)
        real(kind=8),intent(in)     :: dbMean(ndim)
        real(kind=8),intent(in)     :: dbVar(ndim)
        real(kind=8),intent(out)    :: S(ndim)
        real(kind=8),intent(out)    :: D2(ndim)
        real(kind=8),intent(out)    :: Q(ndim)
        real(kind=8),intent(out)    :: DeltaM(ndim)
        real(kind=8),intent(out)    :: DeltaD(ndim)

        integer(kind=4) i, errcode
        real(kind=8) dn

        dn = n
        errcode = 0
        do i = 1, ndim
          D2(i) = c(i,i)*c(i,i)
          Q(i) = 3.0D0*D2(i)
          S(i) = ((Q(i)-D2(i))/dn)-(2*(Q(i)-2*D2(i))/(dn*dn))+           &
     &           ((Q(i)-3*D2(i))/(dn*dn*dn))
          DeltaM(i) = (a(i)-dbMean(i))/sqrt(c(i,i)/dn)
          DeltaD(i) = (c(i,i)-dbVar(i))/sqrt(S(i))
          if (abs(DeltaM(i))>3.0D0 .OR. abs(DeltaD(i))>3.0D0) then
            errcode = 1
          end if
        end do

        dGaussianMVCheckResults = errcode
      END FUNCTION

!
!     Check that single precision GaussianMV sample characteristics agree with theory
!
      FUNCTION sGaussianMVCheckResults( ndim, n, a, c, dbMean, dbVar, S, &
     &      D2, Q, DeltaM, DeltaD )
        integer,intent(in)          :: ndim
        integer,intent(in)          :: n
        real(kind=4),intent(in)     :: a(ndim)
        real(kind=4),intent(in)     :: c(ndim,ndim)
        real(kind=4),intent(in)     :: dbMean(ndim)
        real(kind=4),intent(in)     :: dbVar(ndim)
        real(kind=4),intent(out)    :: S(ndim)
        real(kind=4),intent(out)    :: D2(ndim)
        real(kind=4),intent(out)    :: Q(ndim)
        real(kind=4),intent(out)    :: DeltaM(ndim)
        real(kind=4),intent(out)    :: DeltaD(ndim)

        integer(kind=4) i, errcode
        real(kind=4) dn

        dn = n
        errcode = 0
        do i = 1, ndim
          D2(i) = c(i,i)*c(i,i)
          Q(i) = 3.0*D2(i)
          S(i) = ((Q(i)-D2(i))/dn)-(2*(Q(i)-2*D2(i))/(dn*dn))+           &
     &           ((Q(i)-3*D2(i))/(dn*dn*dn))
          DeltaM(i) = (a(i)-dbMean(i))/sqrt(c(i,i)/dn)
          DeltaD(i) = (c(i,i)-dbVar(i))/sqrt(S(i))
          if (abs(DeltaM(i))>3.0 .OR. abs(DeltaD(i))>3.0) then
            errcode = 1
          end if
        end do

        sGaussianMVCheckResults = errcode
      END FUNCTION
