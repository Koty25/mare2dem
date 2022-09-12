!===============================================================================
! Copyright 2005-2020 Intel Corporation.
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

!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Extended Eigensolvers
!             FORTRAN 77 example
!
!*******************************************************************************
!
! Example program for using Intel MKL Extended Eigensolvers (banded format).
! 
! The following routines are used in the example:
!          SGEMM  SFEAST_SBEV  SFEAST_SBGV  FEASTINIT.
!
! Consider the matrix A 
!
!                 |  5   2   1   1   0   0   0   0   0   0   0  |
!                 |  2   6   3   1   1   0   0   0   0   0   0  |
!                 |  1   3   6   3   1   1   0   0   0   0   0  |
!                 |  1   1   3   6   3   1   1   0   0   0   0  |
!                 |  0   1   1   3   6   3   1   1   0   0   0  |
!    A    =       |  0   0   1   1   3   6   3   1   1   0   0  |,
!                 |  0   0   0   1   1   3   6   3   1   1   0  |
!                 |  0   0   0   0   1   1   3   6   3   1   1  |
!                 |  0   0   0   0   0   1   1   3   6   3   1  | 
!                 |  0   0   0   0   0   0   1   1   3   6   2  |
!                 |  0   0   0   0   0   0   0   1   1   2   5  |
!
!
! stored in LAPACK banded format (DOUBLE PRECISION version).
! B is a unit matrix:
!
!                 |  1   0   0   0   0   0   0   0   0   0   0  |
!                 |  0   1   0   0   0   0   0   0   0   0   0  |
!                 |  0   0   1   0   0   0   0   0   0   0   0  |
!                 |  0   0   0   1   0   0   0   0   0   0   0  |
!                 |  0   0   0   0   1   0   0   0   0   0   0  |
!    B    =       |  0   0   0   0   0   1   0   0   0   0   0  |.
!                 |  0   0   0   0   0   0   1   0   0   0   0  |
!                 |  0   0   0   0   0   0   0   1   0   0   0  |
!                 |  0   0   0   0   0   0   0   0   1   0   0  | 
!                 |  0   0   0   0   0   0   0   0   0   1   0  |
!                 |  0   0   0   0   0   0   0   0   0   0   1  |
!
!  In what follows the symbol ' represents a transpose operation.
!
!  The test performs the following operations :
!
!       1. The code calls  FEASTINIT  to define the default values for the input
!          FEAST parameters.
!
!       2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using 
!          SFEAST_SBEV.
!
!       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
!           are the expected eigenvalues  and E(i) are eigenvalues computed 
!           with the help of SFEAST_SBEV().
!
!       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
!          where X is the matrix of eigenvectors computed with the help of 
!          SFEAST_SBEV. SGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!       5. The  code solves  the generalized eigenvalue problem  Ax=eBx using
!          SFEAST_SBGV.
!
!       6. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
!           are the expected eigenvalues  and E(i) are eigenvalues computed 
!           with the help of SFEAST_SBGV().
!
!       7. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
!          where X is the matrix of eigenvectors computed with the help of 
!          SFEAST_SBEV. SGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!  
!*******************************************************************************
      program sexample_banded
      implicit none
!!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      character*1 UPLO
      parameter  (UPLO='F')
      integer     n, bandA, lda, ldb, bandB
      parameter  (n=11, bandA=3, lda=7, ldb=3, bandB=1  )
      real        A(lda,n), B(ldb,n)
!!!!!!!!!!!!!!!!! Declaration of FEAST variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
      integer     fpm(128)
      real        Emin,Emax
      real        epsout
      integer     loop
      integer     L
      parameter   (L=8)      
      integer     M0,M,info
      real        E(n)
      real        X(n,n)
      real        res(n)
!!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!!!!!!!
      real        Eig(n)
      real        R(n)
      real        Y(n,n)
      integer     i,j
      integer     ldx, ldy
      real        one, zero, smax, eigabs 
!!!!!!!!!!!!!!!!! Exact eigenvalues in range (3.0, 7.0) !!!!!!!!!!!!!!!!!!!!!!
      Eig=0.0
      Eig(1)=3.17157287525381
      Eig(2)=4.0
      Eig(3)=4.0
      Eig(4)=4.1292484841890931
      Eig(5)=4.4066499006731521
      Eig(6)=6.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Initialize input matrices A and B !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      A=0.0
      B=0.0
      
      do i=3,n
        A(1,i)=1.0
        A(2,i)=1.0
        if(i.gt.3)then
          A(6,i-3)=1.0
          A(7,i-3)=1.0
        endif
      enddo
      A(2,3)=1.0
      A(6,9)=1.0

      do i=3,n-1
        A(3,i)=3.0
        A(5,i-1)=3.0
      enddo
      A(3,2)=2.0
      A(3,11)=2.0
      A(5,1)=2.0
      A(5,10)=2.0

      do i=2,n-1
        A(4,i)=6.0
      enddo
      A(4,1)=5.0
      A(4,11)=5.0
      A(5,11)=0.0
      A(6,11)=0.0
      A(1,3)=0.0

      do i=1,n
        B(2,i)=1.0
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!   Print matrix dimension          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print *,' Banded matrix size',n

!!! Search interval [Emin,Emax] including M eigenpairs !!!!!!!!!!!!!!!!!!!!!!!!!
      Emin=3.0
      Emax=7.0
      M0=L
      M=M0
      print *,'Search interval ', Emin,' ', Emax
      ldx=n
      ldy=n
      one=1.0 
      zero=0.0
!
!        Task 1. Call FEASTINIT to define the default values for the input
!        FEAST parameters.
!
      call feastinit(fpm)
      fpm(1)=1
      print *, ' Testing dfeast_sbev '
!
!         Task 2. Solve the standard eigenvalue problem Ax=ex.
! 
      call sfeast_sbev(UPLO,N, bandA,A,lda,fpm,epsout,loop,
     &       Emin,Emax,M0,E,X,M,res,info)
      print  *,' FEAST OUTPUT INFO ',info
      if(info.ne.0) stop 1
!
!         Task 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!         are the expected eigenvalues  and E(i) are eigenvalues computed
!         with the help of SFEAST_SBEV().
!
      print *, 'Number of eigenvalues found ', M
      print *, ' Computed    |    Expected  '
      print *, ' Eigenvalues |    Eigenvalues '
      eigabs=zero
      do i=1,M
         R(i)=abs(E(i)-Eig(i))
         eigabs=max(eigabs, R(i))
         print *, E(i), Eig(i)
      enddo
      print *, ' Max value of '
      print *, ' | computed eigenvalue -expected eigenvalues | ', eigabs
!
!         Task 4.  Compute the maximum absolute value of the matrix
!         Y=(X')*X-I  where X is the matrix of eigenvectors computed with
!         the help of SFEAST_SBEV.
!         Call SGEMM (BLAS Level 3 Routine) to compute (X')*X.
!
      Y=zero
      call sgemm('T','N',M,M,n,one,X,ldx,X,ldx,zero,Y,ldy)
!
!          Compute Y=Y-I
!
       do i=1,M
        Y(i,i)=Y(i,i)-one
      enddo
      print *,'*************************************************'
      print *,'************** REPORT ***************************'
      print *,'*************************************************'
      print *,'# Search interval [Emin,Emax]',Emin,Emax
      print *,'# mode found/subspace',M,M0
      print *,'# iterations',loop
      print *,'Relative error on the trace',epsout
      print *,'Eigenvalues/Residuals'
      do i=1,M
            print *,i,E(i),' | ', res(i)
      enddo
      smax=zero
      do i=1,M
        do j=1,M
          smax=max(smax, abs(Y(i, j)))
        enddo
      enddo
      print *, '  Max absolute value of the matrix '
      print *, ' (transposed of X)*X-I ', smax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! GENERALIZED EIGENVALUE PROBLEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      M0=L
      M=M0
      Y=zero
      E=zero
      X=zero
      res=zero
      epsout=zero
      print *, ' Testing dfeast_sbgv '
!
!         Task 5. Solve the generalized eigenvalue problem Ax=eBx.
!
      call sfeast_sbgv(UPLO,N,bandA, A, lda, bandB, B, ldB,fpm,
     &           epsout,loop,Emin,Emax,M0,E,X,M,res,info)

      print  *,' FEAST OUTPUT INFO ',info
      if(info.ne.0) stop 1
!
!         Task 6. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!         are the expected eigenvalues  and E(i) are eigenvalues computed
!         with the help of SFEAST_SBGV().
!
      print *, 'Number of eigenvalues found ', M
      print *, ' Computed    |    Expected  '
      print *, ' Eigenvalues |    Eigenvalues '
      eigabs=zero
      do i=1,M
         R(i)=abs(E(i)-Eig(i))
         eigabs=max(eigabs, R(i))
         print *, E(i), Eig(i)
      enddo
      print *, ' Max | computed eigenvalue -expected eigenvalues | ',
     & eigabs
!
!         Task 7.  Compute the maximum absolute value of the matrix
!         Y=(X')*X-I  where X is the matrix of eigenvectors computed with
!         the help of SFEAST_SBEV.
!         Call SGEMM (BLAS Level 3 Routine) to compute (X')*X.
!
      Y=zero
      call sgemm('T','N',M,M,n,one,X,ldx,X,ldx,zero,Y,ldy)
!
!          Compute Y=Y-I
!
      do i=1,M
        Y(i,i)=Y(i,i)-one
      enddo
      print *,'*************************************************'
      print *,'************** REPORT ***************************'
      print *,'*************************************************'
      print *,'# Search interval [Emin,Emax]',Emin,Emax
      print *,'# mode found/subspace',M,M0
      print *,'# iterations',loop
      print *,'Relative error on the Trace',epsout
      print *,'Eigenvalues/Residuals'
      do i=1,M
         print *,i,E(i),' | ', res(i)
      enddo
      smax=zero
      do i=1,M
        do j=1,M
          smax=max(smax, abs(Y(i, j)))
        enddo
      enddo
      print *, '  Max absolute value of the matrix '
      print *, ' (transposed of X)*X-I ', smax
      stop 0 
      end program sexample_banded
