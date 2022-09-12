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
! Example program for using Intel MKL Extended Eigensolvers (sparse format).
! 
! The following routines are used in the example:
!          DGEMM  DFEAST_SCSREV  DFEAST_SCSRGV  FEASTINIT.
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
! stored as a compessed sparse row  matrix (DOUBLE PRECISION).
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
!
!  In what follows the symbol ' represents a transpose operation.
!
!  The test performs the following operations :
!
!
!       1. The code calls  FEASTINIT  to define the default values for the input
!          FEAST parameters.
!
!       2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using 
!          DFEAST_SCSREV.
!
!       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
!           are the expected eigenvalues  and E(i) are eigenvalues computed 
!           with the help of DFEAST_SCSREV().
!
!       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
!          where X is the matrix of eigenvectors computed with the help of 
!          DFEAST_SCSREV. DGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!       5. The  code solves  the generalized eigenvalue problem  Ax=eBx using
!          DFEAST_SCSRGV.
!
!       6. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
!           are the expected eigenvalues  and E(i) are eigenvalues computed 
!           with the help of DFEAST_SCSRGV().
!
!       7. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
!          where X is the matrix of eigenvectors computed with the help of 
!          DFEAST_SCSREV. DGEMM (BLAS Level 3 Routine) is called  to compute (X')*X. 
!
!*******************************************************************************
      program  dexample_sparse
      implicit none
!!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      character*1 UPLO
      parameter   (UPLO='F')
      integer     n
      parameter   (n=11)
      real*8      val(65), valb(n)
      integer     rows(n+1), cols(65), rowsb(n+1), colsb(n)
!!!!!!!!!!!!!!!!! Declaration of FEAST variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
      integer     fpm(128)
      real*8      Emin,Emax
      real*8      epsout
      integer     loop
      integer     L
      parameter   (L=8)      
      integer     M0,M,info
      real*8      E(n)
      real*8      X(n,n)
      real*8      res(n)
!!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!!!!!!!
      real*8        Eig(n)
      real*8        R(n)
      real*8        Y(n,n)
      integer       i,j
      integer       ldx, ldy
      real*8        one, zero, smax, eigabs
!!!!!!!!!!!!!!!!! Exact eigenvalues in range (3.0, 7.0) !!!!!!!!!!!!!!!!!!!!!!
      Eig=0.0d0
      Eig(1)=3.17157287525381d0
      Eig(2)=4.0d0
      Eig(3)=4.0d0
      Eig(4)=4.1292484841890931d0
      Eig(5)=4.4066499006731521d0
      Eig(6)=6.0d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Initialize input matrices A and B !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1, 65
         val(i)=1.0d0
      enddo
      do i=1,6
        val(7*(i+1)-2)=6.0d0
        val(7*(i+1)-1)=3.0d0
        val(7*(i+1)+4)=3.0d0
      enddo
      val(1)=5.0d0
      val(2)=2.0d0
      val(5)=2.0d0
      val(6)=6.0d0
      val(7)=3.0d0
      val(11)=3.0d0
      val(54)=6.0d0
      val(55)=3.0d0
      val(59)=3.0d0
      val(60)=6.0d0
      val(61)=2.0d0
      val(64)=2.0d0
      val(65)=5.0d0

      rows(1)=1
      rows(2)=5
      rows(3)=10
      rows(4)=16
      rows(5)=23
      rows(6)=30
      rows(7)=37
      rows(8)=44
      rows(9)=51
      rows(10)=57
      rows(11)=62
      rows(12)=66

      do i=1,5
        do j=2,8
            cols(7*(i+1)+j)=i-1+j-1
        enddo
      enddo

      do j=1,4
        do i=1,2
            cols(5*i+j-1)=j
            cols(48+5*i+j-1)=j+7
        enddo
        cols(j)=j
        cols(61+j)=j+7
      enddo

      cols(9)=5
      cols(14)=5
      cols(15)=6
      cols(51)=6
      cols(52)=7
      cols(57)=7

      do i=1,n
        valb(i)=1.0d0
        colsb(i)=i
        rowsb(i)=i
      enddo
      rowsb(n+1)=n+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!   Print matrix dimension          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print *,'Sparse matrix size',n

!!! Search interval [Emin,Emax] including M eigenpairs !!!!!!!!!!!!!!!!!!!!!!!!!
      Emin=3.0
      Emax=7.0 
      M0=L
      M=M0
      print *,'Search interval ', Emin,' ', Emax
      ldx=n
      ldy=n
      one=1.0d0
      zero=0.0d0

!
!        Task 1. Call FEASTINIT to define the default values for the input
!        FEAST parameters.
!
      call feastinit(fpm)
      fpm(1)=1
      print *, ' Testing dfeast_scsrev '
!
!         Task 2. Solve the standard eigenvalue problem Ax=ex.
! 
      call dfeast_scsrev(UPLO,N,val,rows,cols,fpm,epsout,loop,
     $   Emin,Emax,M0,E,X,M,res,info)
      print  *,' FEAST OUTPUT INFO ',info
      if(info.ne.0) stop 1
!
!         Task 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!         are the expected eigenvalues  and E(i) are eigenvalues computed
!         with the help of SFEAST_SCSREV().
!
      print *, 'Number of eigenvalues found ', M
      print *, ' Computed    |    Expected  '
      print *, ' Eigenvalues |    Eigenvalues '
      eigabs=0.0d0
      do i=1,M
         R(i)=dabs(E(i)-Eig(i))
         eigabs=max(eigabs, R(i))
         print *, E(i), Eig(i)
      enddo
      print *, ' Max value of '
      print *, ' | computed eigenvalue -expected eigenvalues | ', eigabs
!
!         Task 4.  Compute the maximum absolute value of the matrix
!         Y=(X')*X-I  where X is the matrix of eigenvectors computed with
!         the help of DFEAST_SCSREV.
!         Call DGEMM (BLAS Level 3 Routine) to compute (X')*X.
!
      Y=zero
      call dgemm('T','N',M,M,n,one,X,ldx,X,ldx,zero,Y,ldy)
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
      print *, ' Testing dfeast_scsrgv'
!
!         Task 5. Solve the generalized eigenvalue problem Ax=eBx.
!
      call dfeast_scsrgv(UPLO,N,val,rows,cols,valb,rowsb,colsb,
     $ fpm,epsout,loop,Emin,Emax,M0,E,X,M,res,info)
      print  *,' FEAST OUTPUT INFO ',info
      if(info.ne.0) stop 1
!
!         Task 6. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!         are the expected eigenvalues  and E(i) are eigenvalues computed
!         with the help of DFEAST_SCSRGV().
!
      print *, 'Number of eigenvalues found ', M
      print *, ' Computed    |    Expected  '
      print *, ' Eigenvalues |    Eigenvalues '
      eigabs=zero
      do i=1,M
         R(i)=dabs(E(i)-Eig(i))
         eigabs=max(eigabs, R(i))
         print *, E(i), Eig(i)
      enddo
      print *, ' Max | computed eigenvalue -expected eigenvalues | ',
     & eigabs
!
!         Task 7.  Compute the maximum absolute value of the matrix
!         Y=(X')*X-I  where X is the matrix of eigenvectors computed with
!         the help of DFEAST_SCSREV.
!         Call DGEMM (BLAS Level 3 Routine) to compute (X')*X.
!
      Y=zero
      call dgemm('T','N',M,M,n,one,X,ldx,X,ldx,zero,Y,ldy)
!
!          Compute Y=Y-I
!
      do i=1,M
        Y(i,i)=Y(i,i) - one 
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
      stop 0
      end program dexample_sparse
