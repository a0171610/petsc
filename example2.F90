program petsc_mat
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
  PetscErrorCode ierr
  integer, parameter :: n = 3
  real(8) :: x(3), A(3, 3)
  Vec :: xt
  Mat :: At

  x(1) = 2.0d0
  x(2) = 12.0d0
  x(3) = 13.0d0

  A(1, 1) = -1.0d0
  A(2, 3) = 12.0d0

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call transform_to_petsc_vector(x, xt, n, ierr)
  call transform_to_petsc_matrix(A, At, n, ierr)
   
end program petsc_mat

subroutine transform_to_petsc_vector(x, xt, n, ierr)
  use petscksp
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: x(n)
  integer :: idx(n), i
  PetscScalar :: v(n)
  Vec :: xt
  PetscErrorCode ierr

  do i = 1, n
    idx(i) = i - 1
    v(i) = x(i)
  enddo

  call VecCreateSeq(PETSC_COMM_SELF, n, xt, ierr)
  call VecSetValues(xt, n, idx, v, INSERT_VALUES,ierr)

  call VecView(xt, PETSC_VIEWER_STDOUT_SELF, ierr)

end subroutine transform_to_petsc_vector

subroutine transform_to_petsc_matrix(A, At, n, ierr)
  use petscksp
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: A(n, n)
  Mat :: At
  integer :: i, j, idx(n)
  PetscErrorCode ierr


  do i = 1, n
    idx(i) = i - 1
  enddo

  call MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, n, PETSC_NULL_INTEGER, At, ierr)
  do i = 1, n
    call MatSetValues(At, 1, i - 1, n, idx, A(i, :), INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(At, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(At, MAT_FINAL_ASSEMBLY, ierr)

  call MatView(At, PETSC_VIEWER_STDOUT_SELF, ierr)

end subroutine transform_to_petsc_matrix
