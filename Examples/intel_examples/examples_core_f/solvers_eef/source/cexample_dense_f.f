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
! Example program for using Intel MKL Extended Eigensolvers (dense format).
! 
! The following routines are used in the example:
!          CGEMM  CFEAST_HEEV  SFEAST_HEGV  FEASTINIT.
!
! Consider the matrix A 
!
!                 |  10   1+2i  0    0    0    0    0    0    0    0    |
!                 |  1-2i  9   2+3i  0    0    0    0    0    0    0    |
!                 |  0    2-3i  8   3+4i  0    0    0    0    0    0    |
!                 |  0    0    3-4i  7   4+5i  0    0    0    0    0    |
!                 |  0    0    0    4-5i  6   5+6i  0    0    0    0    |
!    A    =       |  0    0    0    0    5-6i  5   6+7i  0    0    0    |,
!                 |  0    0    0    0    0    6-7i  4   7+8i  0    0    |
!                 |  0    0    0    0    0    0    7-8i  3   8+9i  0    |
!                 |  0    0    0    0    0    0    0    8-9i  2   9+10i | 
!                 |  0    0    0    0    0    0    0    0    9-10i  1   |
!
!
! stored as dense matrix (single complex precision version).
! 
! B is a unit matrix:
!
!                 |  1   0   0   0   0   0   0   0   0   0  |
!                 |  0   1   0   0   0   0   0   0   0   0  |
!                 |  0   0   1   0   0   0   0   0   0   0  |
!                 |  0   0   0   1   0   0   0   0   0   0  |
!                 |  0   0   0   0   1   0   0   0   0   0  |
!    B    =       |  0   0   0   0   0   1   0   0   0   0  |.
!                 |  0   0   0   0   0   0   1   0   0   0  |
!                 |  0   0   0   0   0   0   0   1   0   0  |
!                 |  0   0   0   0   0   0   0   0   1   0  | 
!                 |  0   0   0   0   0   0   0   0   0   1  |
!
!  In what follows the symbol ' represents a conjugate transpose operation.
!
!  The test performs the following operations :
!  
!       1. The code calls  FEASTINIT  to define the default values for the input
!          FEAST parameters.
!
!       2. The  code  solves  the  standard eigenvalue  problem  Ax=ex   using 
!          CFEAST_HEEV.
!
!       3. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
!           are the expected eigenvalues  and E(i) are eigenvalues computed 
!           with the help of CFEAST_HEEV().
!
!       4. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
!          where X is the matrix of eigenvectors computed with the help of 
!          CFEAST_HEEV. CGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!       5. The  code solves  the generalized eigenvalue problem  Ax=eBx using
!          CFEAST_HEGV.
!
!       6. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
!           are the expected eigenvalues  and E(i) are eigenvalues computed 
!           with the help of CFEAST_HEGV().
!
!       7. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
!          where X is the matrix of eigenvectors computed with the help of 
!          CFEAST_HEGV. CGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
! 
!*******************************************************************************
      program cexample_dense
      implicit none
!!!!!!!!!!!!!!!!! Matrix declaration variable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      character*1 UPLO
      parameter   (UPLO='F')
      integer     n
      parameter   (n=10)
      complex     A(n,n), B(n,n)
!!!!!!!!!!!!!!!!! Feast declaration variable !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
      integer     feastparam(128)
      real        Emin,Emax
      real        epsout
      integer     loop
      integer     L
      parameter   (L=6)
      real        E(n)
      complex     X(n,n)
      real        res(n)
!!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!!!!!!!
      real        Eig(n)
      real        R(n)
      complex     Y(n,n)
      integer     i,j
      integer     ldx, ldy
      complex     one, zero
      real        smax, eigabs
      integer     M0,M,info

!!!!!!!!!!!!!!!!! Exact eigenvalues in range (2.0, 12.0)  !!!!!!!!!!!!!!!!!!!!!!
      Eig=0.0
      Eig(1)=2.231051
      Eig(2)=6.058517
      Eig(3)=9.109751
      Eig(4)=11.703148

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Initialize input matrices A and B !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      one=cmplx(1.0,0.0) 
      zero=cmplx(0.0, 0.0)
      A=zero
      B=zero
      do i=1,n-1       
        A(i,i)=dcmplx(dble(n+1-i),0.0)
        A(i,i+1)=dcmplx(dble(i),dble(i+1))  
        A(i+1,i)=dcmplx(dble(i),dble(-(i+1)))
      enddo
      A(10,10)=one



      do i=1,n
        B(i,i)=one
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print *,'Dense matrix size',n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!   Print matrix dimension          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Search interval [Emin,Emax] including M eigenpairs !!!!!!!!!!!!!!!!!!!!!!!!!
      Emin=2.0
      Emax=12.0 
      M0=L
      M=M0
      print *,'Search interval ', Emin,' ', Emax
      ldx=n
      ldy=n
!
!        Task 1. Call FEASTINIT to define the default values for the input
!        FEAST parameters.
!
      call feastinit(feastparam)
      feastparam(1)=1
      print *, ' Testing cfeast_heev '
!
!         Task 2. Solve the standard eigenvalue problem Ax=ex.
!
      call cfeast_heev(UPLO,n, A,n, feastparam, epsout, loop,Emin,Emax,
     $ M0,E,X,M,res,info)
      print  *,' FEAST OUTPUT INFO ',info
      if(info.ne.0) stop 1
!
!         Task 3. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!         are the expected eigenvalues  and E(i) are eigenvalues computed
!         with the help of CFEAST_HEEV().
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
!         the help of CFEAST_HEEV.
!         Call CGEMM (BLAS Level 3 Routine) to compute (X')*X.
!
      Y=zero
      call cgemm('C','N',M,M,n,one,X,ldx,X,ldx,zero,Y,ldy)
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
      print *, '  Max  value of '
      print *, ' (conjugate transposed of X)*X-I ', smax
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
      print *, ' Testing cfeast_hegv '
!
!         Task 5. Solve the generalized eigenvalue problem Ax=eBx.
!
      call cfeast_hegv(UPLO,N,A(1,1),N,B,N,feastparam,epsout,loop,Emin,
     $ Emax,M0,E,X,M,res,info)
      print  *,' FEAST OUTPUT INFO ',info
      if(info.ne.0) stop 1
!
!         Task 6. Compute the residual R(i) = | E(i) - Eig(i) | where Eig(i)
!         are the expected eigenvalues  and E(i) are eigenvalues computed
!         with the help of CFEAST_HEGV().
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
!         the help of CFEAST_HEEV.
!         Call CGEMM (BLAS Level 3 Routine) to compute (X')*X.
!
      Y=zero
      call cgemm('C','N',M,M,n,one,X,ldx,X,ldx,zero,Y,ldy)
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
      print *, '  Max  value of '
      print *, ' (conjugate transposed of X)*X-I ', smax
      stop 0
      end program cexample_dense
