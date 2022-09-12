!===============================================================================
! Copyright 2020 Intel Corporation.
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

module common_blas
contains 

  subroutine dinit_matrix (trans, row, col, ld, mat)
    character*1,intent(in)         :: trans
    integer,intent(in)             :: row, col, ld
    double precision,intent(inout) :: mat(ld,*)

    integer :: i, j
    
    do i = 1, row
       do j = 1, col
          if ((trans.eq.'N').or.(trans.eq.'n')) then
             mat(i,j) = DBLE(i + ld * j) / DBLE(ld * col)
          else
             mat(j,i) = DBLE(j + ld * i) / DBLE(ld * row)
          end if
       end do
    end do
  end subroutine dinit_matrix

  subroutine dcopy_matrix (trans, row, col, ld, src, dest)
    character*1,intent(in)         :: trans
    integer,intent(in)             :: row, col, ld
    double precision,intent(in)    :: src(ld,*)
    double precision,intent(inout) :: dest(ld,*)

    integer :: i, j
    
    do i = 1, row
       do j = 1, col
          dest(i,j) = src(i,j)
       end do
    end do
  end subroutine dcopy_matrix

  function dcheck_matrix (row, col, ld, res, ref)
    integer                        :: dcheck_matrix
    integer,intent(in)             :: row, col, ld
    double precision,intent(in)    :: ref(ld,*), res(ld,*)

    integer          :: i, j
    double precision :: diff
    
    dcheck_matrix = 0
    do i = 1, row
       do j = 1, col
          diff = abs(res(i,j) - ref(i,j))
          if (diff.gt.0.00001) then
             print 100, i, j, ref(i,j), res(i,j), diff 
             dcheck_matrix = 1
          end if
       end do
    end do

100 format(7x, 'Error at index (',i1,',',i1,'), expected = ',f10.6,', computed = ',f10.6,', difference = ',f10.6)
  end function dcheck_matrix

  subroutine dprint_matrix (row, col, ld, mat)
    integer,intent(in)             :: row, col
    double precision,intent(in)    :: mat(ld,col)

    integer :: i

    do i = 1, row
       write( * , '(10G14.6)') mat(i,:)
    end do

  end subroutine dprint_matrix

  subroutine sinit_matrix (trans, row, col, ld, mat)
    character*1,intent(in)         :: trans
    integer,intent(in)             :: row, col, ld
    real,intent(inout) :: mat(ld,*)

    integer :: i, j
    
    do i = 1, row
       do j = 1, col
          if ((trans.eq.'N').or.(trans.eq.'n')) then
             mat(i,j) = REAL(i + ld * j) / REAL(ld * col)
          else
             mat(j,i) = REAL(j + ld * i) / REAL(ld * row)
          end if
       end do
    end do
  end subroutine sinit_matrix

  subroutine scopy_matrix (trans, row, col, ld, src, dest)
    character*1,intent(in)         :: trans
    integer,intent(in)             :: row, col, ld
    real,intent(in)    :: src(ld,*)
    real,intent(inout) :: dest(ld,*)

    integer :: i, j
    
    do i = 1, row
       do j = 1, col
          dest(i,j) = src(i,j)
       end do
    end do
  end subroutine scopy_matrix

  function scheck_matrix (row, col, ld, res, ref)
    integer                        :: scheck_matrix
    integer,intent(in)             :: row, col, ld
    real,intent(in)    :: ref(ld,*), res(ld,*)

    integer          :: i, j
    real :: diff
    
    scheck_matrix = 0
    do i = 1, row
       do j = 1, col
          diff = abs(res(i,j) - ref(i,j))
          if (diff.gt.0.0001) then
             print 100, i, j, ref(i,j), res(i,j), diff 
             scheck_matrix = 1
          end if
       end do
    end do

100 format(7x, 'Error at index (',i1,',',i1,'), expected = ',f10.6,', computed = ',f10.6,', difference = ',f10.6)
  end function scheck_matrix

  subroutine sprint_matrix (row, col, ld, mat)
    integer,intent(in)             :: row, col
    real,intent(in)    :: mat(ld,col)

    integer :: i

    do i = 1, row
       write( * , '(10G14.6)') mat(i,:)
    end do

  end subroutine sprint_matrix

  subroutine cinit_matrix (trans, row, col, ld, mat)
    character*1,intent(in)         :: trans
    integer,intent(in)             :: row, col, ld
    complex,intent(inout) :: mat(ld,*)

    integer :: i, j
    
    do i = 1, row
       do j = 1, col
          if ((trans.eq.'N').or.(trans.eq.'n')) then
             mat(i,j) = CMPLX(DBLE(i + ld * j) / DBLE(ld * col), &
                  & 1 - DBLE(i + ld * j) / DBLE(ld * col))
          else
             mat(j,i) = CMPLX(DBLE(j + ld * i) / DBLE(ld * row), &
                  & 1 - DBLE(j + ld * i) / DBLE(ld * row))
          end if
       end do
    end do
  end subroutine cinit_matrix

  subroutine ccopy_matrix (trans, row, col, ld, src, dest)
    character*1,intent(in)         :: trans
    integer,intent(in)             :: row, col, ld
    complex,intent(in)    :: src(ld,*)
    complex,intent(inout) :: dest(ld,*)

    integer :: i, j
    
    do i = 1, row
       do j = 1, col
          dest(i,j) = src(i,j)
       end do
    end do
  end subroutine ccopy_matrix

  function ccheck_matrix (row, col, ld, res, ref)
    integer                        :: ccheck_matrix
    integer,intent(in)             :: row, col, ld
    complex,intent(in)    :: ref(ld,*), res(ld,*)

    integer          :: i, j
    real :: diff
    
    ccheck_matrix = 0
    do i = 1, row
       do j = 1, col
          diff = abs(res(i,j) - ref(i,j))
          if (diff.gt.0.0001) then
             print 100, i, j, real(ref(i,j)), imag(ref(i,j)), real(res(i,j)), imag(res(i,j)), diff
             ccheck_matrix = 1
          end if
       end do
    end do

100 format(7x, 'Error at index (',i1,',',i1,'), expected = ',f10.6,' + i ',f10.6,', computed = ',f10.6,' + i ',f10.6,', difference = ',f10.6)
  end function ccheck_matrix

  subroutine cprint_matrix (row, col, ld, mat)
    integer,intent(in)             :: row, col
    complex,intent(in)    :: mat(ld,col)

    integer :: i

    do i = 1, row
       write( * , '(10G14.6)') mat(i,:)
    end do

  end subroutine cprint_matrix

  subroutine zinit_matrix (trans, row, col, ld, mat)
    character*1,intent(in)         :: trans
    integer,intent(in)             :: row, col, ld
    complex*16,intent(inout) :: mat(ld,*)

    integer :: i, j
    
    do i = 1, row
       do j = 1, col
          if ((trans.eq.'N').or.(trans.eq.'n')) then
             mat(i,j) = CMPLX(DBLE(i + ld * j) / DBLE(ld * col), &
                  & 1 - DBLE(i + ld * j) / DBLE(ld * col), 16)
          else
             mat(j,i) = CMPLX(DBLE(j + ld * i) / DBLE(ld * row), &
                  & 1 - DBLE(j + ld * i) / DBLE(ld * row), 16)
          end if
       end do
    end do
  end subroutine zinit_matrix

  subroutine zcopy_matrix (trans, row, col, ld, src, dest)
    character*1,intent(in)         :: trans
    integer,intent(in)             :: row, col, ld
    complex*16,intent(in)    :: src(ld,*)
    complex*16,intent(inout) :: dest(ld,*)

    integer :: i, j
    
    do i = 1, row
       do j = 1, col
          dest(i,j) = src(i,j)
       end do
    end do
  end subroutine zcopy_matrix

  function zcheck_matrix (row, col, ld, res, ref)
    integer                        :: zcheck_matrix
    integer,intent(in)             :: row, col, ld
    complex*16,intent(in)    :: ref(ld,*), res(ld,*)

    integer          :: i, j
    double precision :: diff
    
    zcheck_matrix = 0
    do i = 1, row
       do j = 1, col
          diff = abs(res(i,j) - ref(i,j))
          if (diff.gt.0.00001) then
             print 100, i, j, real(ref(i,j)), imag(ref(i,j)), real(res(i,j)), imag(res(i,j)), diff
             zcheck_matrix = 1
          end if
       end do
    end do

100 format(7x, 'Error at index (',i1,',',i1,'), expected = ',f10.6,' + i ',f10.6,', computed = ',f10.6,' + i ',f10.6,', difference = ',f10.6)
  end function zcheck_matrix

  subroutine zprint_matrix (row, col, ld, mat)
    integer,intent(in)             :: row, col
    complex*16,intent(in)    :: mat(ld,col)

    integer :: i

    do i = 1, row
       write( * , '(10G14.6)') mat(i,:)
    end do

  end subroutine zprint_matrix

  subroutine dinit_vector (n, inc, vec)
    integer,intent(in)             :: n, inc
    double precision,intent(inout) :: vec(*)

    integer :: i, idx

    idx = 1

    if (inc.lt.0) idx = 1 + (1 - n) * inc
    
    do i = 1, n
       vec(idx) = DBLE(idx) / DBLE(1 + (n - 1) * abs(inc))
       idx = idx + inc
    end do
  end subroutine dinit_vector

  subroutine dcopy_vector (n, inc, src, dest)
    integer,intent(in)             :: n, inc
    double precision,intent(in)    :: src(*)
    double precision,intent(inout) :: dest(*)

    integer :: i, idx

    idx = 1

    if (inc.lt.0) idx = 1 + (1 - n) * inc
    
    do i = 1, n
          dest(idx) = src(idx)
       idx = idx + inc
    end do
    
  end subroutine dcopy_vector

  function dcheck_vector (n, inc, res, ref)
    integer                        :: dcheck_vector
    integer,intent(in)             :: n, inc
    double precision,intent(in)    :: ref(*), res(*)

    integer :: i, idx
    double precision :: diff

    idx = 1
    
    if (inc.lt.0) idx = 1 + (1 - n) * inc

    dcheck_vector = 0
    
    do i = 1, n
       diff = abs(res(idx) - ref(idx))
       if (diff.gt.0.00001) then
          print 100, i, ref(idx), res(idx), diff 
          dcheck_vector = 1
       end if
       idx = idx + inc
    end do
    
100 format(7x, 'Error at index (',i1,'), expected = ',f10.6,', computed = ',f10.6,', difference = ',f10.6)
  end function dcheck_vector

  
  subroutine sinit_vector (n, inc, vec)
    integer,intent(in)             :: n, inc
    real,intent(inout) :: vec(*)

    integer :: i, idx

    idx = 1

    if (inc.lt.0) idx = 1 + (1 - n) * inc
    
    do i = 1, n
       vec(idx) = REAL(idx) / REAL(1 + (n - 1) * abs(inc))
       idx = idx + inc
    end do
  end subroutine sinit_vector

  subroutine scopy_vector (n, inc, src, dest)
    integer,intent(in)             :: n, inc
    real,intent(in)    :: src(*)
    real,intent(inout) :: dest(*)

    integer :: i, idx

    idx = 1

    if (inc.lt.0) idx = 1 + (1 - n) * inc
    
    do i = 1, n
          dest(idx) = src(idx)
       idx = idx + inc
    end do
    
  end subroutine scopy_vector

  function scheck_vector (n, inc, res, ref)
    integer                        :: scheck_vector
    integer,intent(in)             :: n, inc
    real,intent(in)    :: ref(*), res(*)

    integer :: i, idx
    real :: diff

    idx = 1
    
    if (inc.lt.0) idx = 1 + (1 - n) * inc

    scheck_vector = 0
    
    do i = 1, n
       diff = abs(res(idx) - ref(idx))
       if (diff.gt.0.0001) then
          print 100, i, ref(idx), res(idx), diff 
          scheck_vector = 1
       end if
       idx = idx + inc
    end do
    
100 format(7x, 'Error at index (',i1,'), expected = ',f10.6,', computed = ',f10.6,', difference = ',f10.6)
  end function scheck_vector

  subroutine zinit_vector (n, inc, vec)
    integer,intent(in)             :: n, inc
    complex*16,intent(inout) :: vec(*)

    integer :: i, idx

    idx = 1

    if (inc.lt.0) idx = 1 + (1 - n) * inc
    
    do i = 1, n
       vec(idx) = CMPLX(DBLE(idx) / DBLE(1 + (n - 1) * abs(inc)), &
            & 1 - (DBLE(idx) / DBLE(1 + (n - 1) * abs(inc))), 16)
       idx = idx + inc
    end do
  end subroutine zinit_vector

  subroutine zcopy_vector (n, inc, src, dest)
    integer,intent(in)             :: n, inc
    complex*16,intent(in)    :: src(*)
    complex*16,intent(inout) :: dest(*)

    integer :: i, idx

    idx = 1

    if (inc.lt.0) idx = 1 + (1 - n) * inc
    
    do i = 1, n
          dest(idx) = src(idx)
       idx = idx + inc
    end do
    
  end subroutine zcopy_vector

  function zcheck_vector (n, inc, res, ref)
    integer                        :: zcheck_vector
    integer,intent(in)             :: n, inc
    complex*16,intent(in)    :: ref(*), res(*)

    integer :: i, idx
    double precision :: diff

    idx = 1
    
    if (inc.lt.0) idx = 1 + (1 - n) * inc

    zcheck_vector = 0
    
    do i = 1, n
       diff = abs(res(idx) - ref(idx))
       if (diff.gt.0.00001) then
          print 100, i, real(ref(idx)), imag(ref(idx)), real(res(idx)), imag(res(idx)), diff
          zcheck_vector = 1
       end if
       idx = idx + inc
    end do
    
100 format(7x, 'Error at index (',i1,'), expected = ',f10.6,' + i ',f10.6,', computed = ',f10.6,' + i ',f10.6,', difference = ',f10.6)
  end function zcheck_vector

  subroutine cinit_vector (n, inc, vec)
    integer,intent(in)             :: n, inc
    complex,intent(inout) :: vec(*)

    integer :: i, idx

    idx = 1

    if (inc.lt.0) idx = 1 + (1 - n) * inc
    
    do i = 1, n
       vec(idx) = CMPLX(DBLE(idx) / DBLE(1 + (n - 1) * abs(inc)), &
            & 1 - (DBLE(idx) / DBLE(1 + (n - 1) * abs(inc))))
       idx = idx + inc
    end do
  end subroutine cinit_vector

  subroutine ccopy_vector (n, inc, src, dest)
    integer,intent(in)             :: n, inc
    complex,intent(in)    :: src(*)
    complex,intent(inout) :: dest(*)

    integer :: i, idx

    idx = 1

    if (inc.lt.0) idx = 1 + (1 - n) * inc
    
    do i = 1, n
          dest(idx) = src(idx)
       idx = idx + inc
    end do
    
  end subroutine ccopy_vector

  function ccheck_vector (n, inc, res, ref)
    integer                        :: ccheck_vector
    integer,intent(in)             :: n, inc
    complex,intent(in)    :: ref(*), res(*)

    integer :: i, idx
    real :: diff

    idx = 1
    
    if (inc.lt.0) idx = 1 + (1 - n) * inc

    ccheck_vector = 0
    
    do i = 1, n
       diff = abs(res(idx) - ref(idx))
       if (diff.gt.0.0001) then
          print 100, i, real(ref(idx)), imag(ref(idx)), real(res(idx)), imag(res(idx)), diff
          ccheck_vector = 1
       end if
       idx = idx + inc
    end do
    
100 format(7x, 'Error at index (',i1,'), expected = ',f10.6,' + i ',f10.6,', computed = ',f10.6,' + i ',f10.6,', difference = ',f10.6)
  end function ccheck_vector

  function dcheck_scalar (res, ref)
    integer                        :: dcheck_scalar
    double precision,intent(in)    :: ref, res

    double precision :: diff
    
    dcheck_scalar = 0
    
    diff = abs(res - ref)
    if (diff.gt.0.00001) then
       print 100, ref, res, diff 
       dcheck_scalar = 1
    end if
    
100 format(7x, 'Error expected = ',f10.6,', computed = ',f10.6,', difference = ',f10.6)
  end function dcheck_scalar

  function scheck_scalar (res, ref)
    integer                        :: scheck_scalar
    real,intent(in)    :: ref, res

    real :: diff
    
    scheck_scalar = 0
    
    diff = abs(res - ref)
    if (diff.gt.0.00001) then
       print 100, ref, res, diff 
       scheck_scalar = 1
    end if
    
100 format(7x, 'Error expected = ',f10.6,', computed = ',f10.6,', difference = ',f10.6)
  end function scheck_scalar

  function zcheck_scalar (res, ref)
    integer                        :: zcheck_scalar
    complex*16,intent(in)    :: ref, res

    double precision :: diff
    
    zcheck_scalar = 0
    
    diff = abs(res - ref)
    if (diff.gt.0.00001) then
       print 100, real(ref), imag(ref), real(res), imag(res), diff 
       zcheck_scalar = 1
    end if
    
100 format(7x, 'Error expected = ',f10.6,' + i ',f10.6,', computed = ',f10.6,' + i ',f10.6,', difference = ',f10.6)
  end function zcheck_scalar

  function ccheck_scalar (res, ref)
    integer                        :: ccheck_scalar
    complex,intent(in)    :: ref, res

    real :: diff
    
    ccheck_scalar = 0
    
    diff = abs(res - ref)
    if (diff.gt.0.00001) then
       print 100, real(ref), imag(ref), real(res), imag(res), diff 
       ccheck_scalar = 1
    end if
    
100 format(7x, 'Error expected = ',f10.6,' + i ',f10.6,', computed = ',f10.6,' + i ',f10.6,', difference = ',f10.6)
  end function ccheck_scalar

end module common_blas
