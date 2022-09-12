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
! Example program for finding largest eigenvalues of standard and generalized
! eigenvalue problem using Intel MKL Extended Eigensolvers (sparse format).
! 
! The following routines are used in the example:
!          SGEMM  MKL_SPARSE_S_EV  MKL_SPARSE_EE_INIT MKL_SPARSE_S_CREATE_CSR
!          MKL_SPARSE_S_GV.
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
! and unit matrix B:
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
!
!  The test performs the following operations :
!
!
!       1. The code creates sparse matrix structure using Sparse BLAS routine
!          MKL_SPARSE_S_CREATE_CSR.
!
!       2. The code calls  MKL_SPARSE_EE_INIT  to define the default values for the input
!          MKL_SPARSE_S_EV parameters.
!
!       3. The  code  finds largest eigenvalues of standard eigenvalue  problem  Ax=ex using 
!          MKL_SPARSE_S_EV.
!
!       4. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
!           are the expected eigenvalues  and E(i) are eigenvalues computed 
!           with the help of MKL_SPARSE_S_EV().
!
!       5. The code computes the maximum absolute value of the matrix  Y=(X')*X-I  
!          where X is the matrix of eigenvectors computed with the help of 
!          MKL_SPARSE_S_EV. SGEMM (BLAS Level 3 Routine) is called  to compute (X')*X.
!
!       6. The  code  finds largest eigenvalues of standard eigenvalue  problem  Ax=eBx using 
!          MKL_SPARSE_S_GV.
!
!       7. The code computes the residual R(i) = | E(i) - Eig(i) |  where Eig(i) 
!           are the expected eigenvalues  and E(i) are eigenvalues computed 
!           with the help of MKL_SPARSE_S_GV().
!
!       8. The  code  finds largest singular values of A using MKL_SPARSE_S_SVD. 
!
!       9. The code computes the residual R(i) = | sigma(i) - Eig(i) |  where Eig(i) 
!           are the expected singula values  and sigma(i) are singular values computed 
!           with the help of MKL_SPARSE_S_SVD().
!
!
!
!*******************************************************************************
PROGRAM  sexample_extremal_ev
      USE MKL_SPBLAS
      USE MKL_SOLVERS_EE
      implicit none
!!!!!!!!!!!!!!!!! Matrix declaration variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer      n
      parameter   (n=11)
      real        val(65), valb(n)
      integer     rows(n+1), cols(65), rowsb(n+1), colsb(n)
      !   Matrix descriptor
      TYPE(MATRIX_DESCR) descrA
      TYPE(MATRIX_DESCR) descrB
      !   CSR matrix structure
      TYPE(SPARSE_MATRIX_T) csrA
      TYPE(SPARSE_MATRIX_T) csrB
!!!!!!!!!!!!!!!!! Declaration of MKL_SPARSE_S_EV variables !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! E - eigenvalues, X - eigenvectors, res - residual !!!!!!!!!!!!
!!!!!!!!!!!!!!!!! XL - left singular vectors, XR - right singular vectors !!!!!!
      character*1 WHICH
      parameter   (WHICH='L') !L - for largest eigenvalues to find
      character*1 WHICHV      
      parameter   (WHICHV='R') !R - for right singular values to find
      integer     pm(128)
      integer     K0,K,info
      parameter   (K0 = 5) !Required number of eigenvalues/singular values to find
      real        E(n)
      real        sigma(n)
      real        X(n,n), XL(n,n), XR(n,n)
      real        res(n)
!!!!!!!!!!!!!!!!! Declaration of local variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Eig - exact eigenvalues, R=|E-Eig|, Y=(X')*X-I !!!!!!!!!!!!!!!
      real        Eig(n)
      real        R(n)
      real        Y(n,n)
      integer     i,j
      integer     ldx, ldy
      real        one, zero, smax
!!!!!!!!!!!!!!!!! Exact eigenvalues  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Eig=0.0
      Eig(1)=0.52228228746139
      Eig(2)=1.80384757729340
      Eig(3)=3.17157287525381
      Eig(4)=4.0
      Eig(5)=4.0
      Eig(6)=4.1292484841890931
      Eig(7)=4.4066499006731521
      Eig(8)=6.0
      Eig(9)=8.82842712474622
      Eig(10)=12.1961524227068
      Eig(11)=14.9418193276764

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Initialize input matrix A !!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1, 65
         val(i)=1.0
      enddo
      do i=1,6
        val(7*(i+1)-2)=6.0
        val(7*(i+1)-1)=3.0
        val(7*(i+1)+4)=3.0
      enddo
      val(1)=5.0
      val(2)=2.0
      val(5)=2.0
      val(6)=6.0
      val(7)=3.0
      val(11)=3.0
      val(54)=6.0
      val(55)=3.0
      val(59)=3.0
      val(60)=6.0
      val(61)=2.0
      val(64)=2.0
      val(65)=5.0

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Initialize input matrix B !!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,n
        valb(i)=1.0
        colsb(i)=i
        rowsb(i)=i
      enddo
      rowsb(n+1)=n+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!   Print matrix dimension          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print *,'Sparse matrix size',n
      ldx=n
      ldy=n
      one=1.0
      zero=0.0
!
!        Task 1. Call MKL_SPARSE_C_CREATE_CSR to create matrix handle
!      

      info = mkl_sparse_s_create_csr(csrA,SPARSE_INDEX_BASE_ONE,N,N,rows,rows(2),cols,val)
      info = mkl_sparse_s_create_csr(csrB,SPARSE_INDEX_BASE_ONE,N,N,rowsb,rowsb(2),colsb,valb)

!         Create matrix descriptor
      descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL
      descrB % TYPE = SPARSE_MATRIX_TYPE_GENERAL
!
!        Task 2. Call MKL_SPARSE_EE_INIT to define the default values for the input
!        parameters.
!
      info = mkl_sparse_ee_init(pm)
      pm(2) = 4 !Setting tolerance
      pm(8) = 1 ! Use absolute stopping criteria
      
      print *, ' Testing mkl_sparse_s_ev '
!
!         Task 3. Solve the standard eigenvalue problem Ax=ex.
! 
      info = mkl_sparse_s_ev(WHICH,pm,csrA,descrA,k0,k,E,X,res)
      print  *,' OUTPUT INFO ',info
      if(info.ne.0) stop 1
!
!         Task 4. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!         are the expected eigenvalues  and E(i) are eigenvalues computed
!         with the help of mkl_sparse_s_ev().
!
      print *, 'Number of eigenvalues found ', k
      print *, ' Computed    |    Expected  '
      print *, ' Eigenvalues |    Eigenvalues '
      do i=1,K
         R(i)=abs(E(i)-Eig((N-K)+i))
         print *, E(i), Eig((N-K)+i)
      enddo
!
!         Task 5.  Compute the maximum absolute value of the matrix
!         Y=(X')*X-I  where X is the matrix of eigenvectors computed with
!         the help of mkl_sparse_s_ev
!         Call SGEMM (BLAS Level 3 Routine) to compute (X')*X.
!
      Y=zero
      call sgemm('T','N',K,K,n,one,X,ldx,X,ldx,zero,Y,ldy)
!
!          Compute Y=Y-I
!
       do i=1,K
        Y(i,i)=Y(i,i)-one
      enddo
      print *,'*************************************************'
      print *,'************** REPORT ***************************'
      print *,'*************************************************'
      print *,'# Number of eigenvalues requested/found',K0,K
      print *,'Eigenvalues/Residuals'
      do i=1,K
        print *,i,E(i),' | ', res(i)
      enddo
      smax=zero
      do i=1,K
        do j=1,K
           smax=max(smax, abs(Y(i, j)))
        enddo
      enddo
      print *, '  Max absolute value of the matrix '
      print *, ' (transposed of X)*X-I ', smax
      print *, ' Testing mkl_sparse_s_gv '
!
!         Task 6. Solve the standard eigenvalue problem Ax=eBx.
! 
      info = mkl_sparse_s_gv(WHICH,pm,csrA,descrA,csrB,descrB,k0,k,E,X,res)
      print  *,' OUTPUT INFO ',info
      if(info.ne.0) stop 1
!
!         Task 7. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!         are the expected eigenvalues  and E(i) are eigenvalues computed
!         with the help of mkl_sparse_s_ev().
!
      print *, 'Number of eigenvalues found ', k
      print *, ' Computed    |    Expected  '
      print *, ' Eigenvalues |    Eigenvalues '
      do i=1,K
         R(i)=abs(E(i)-Eig((N-K)+i))
         print *, E(i), Eig((N-K)+i)
      enddo

      print *,'*************************************************'
      print *,'************** REPORT ***************************'
      print *,'*************************************************'
      print *,'# Number of eigenvalues requested/found',K0,K
      print *,'Eigenvalues/Residuals'
      do i=1,K
        print *,i,E(i),' | ', res(i)
      enddo

      print *, ' Testing mkl_sparse_s_svd '
!
!         Task 8. Find singular values of A.
! 
      info = mkl_sparse_s_svd(WHICH,WHICHV,pm,csrA,descrA,k0,k,sigma,XL,XR,res)
      print  *,' OUTPUT INFO ',info
      if(info.ne.0) stop 1
!
!         Task 9. Compute the residual R(i) = | E(i) - Eig(i) |  where Eig(i)
!         are the expected singular values  and E(i) are singular values computed
!         with the help of mkl_sparse_s_ev().
!
      print *, 'Number of singular values found ', k
      print *, ' Computed        |    Expected  '
      print *, ' Singular values |    Singular values '
      do i=1,K
         R(i)=abs(sigma(i)-Eig((N-K)+i))
         print *, sigma(i), Eig((N-K)+i)
      enddo

      print *,'*************************************************'
      print *,'************** REPORT ***************************'
      print *,'*************************************************'
      print *,'# Number of singular values requested/found',K0,K
      print *,'singular values/Residuals'
      do i=1,K
        print *,i,sigma(i),' | ', res(i)
      enddo

      !   Release internal representation of CSR matrix
      info = MKL_SPARSE_DESTROY(csrA)
      info = MKL_SPARSE_DESTROY(csrB)

end program sexample_extremal_ev
