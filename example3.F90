program petsc_test
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
   
      integer, parameter :: n = 4
      integer :: ierr, idx(n)
      Vec :: v
      PetscScalar, pointer :: p_v(:)
      PetscScalar :: val(n), x
   
      val(:) = 1d0; val(3) = 4d0;
      idx = (/0,1,2,3/) !0 start indexに注意
   
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      call VecCreateSeq(PETSC_COMM_SELF, n, v, ierr)
   
      call VecSetValues(v, n, idx, val, INSERT_VALUES, ierr)
      x = 3d0 !手間だが一旦格納する
      call VecSetValue(v, 1, x, INSERT_VALUES,ierr) !idx:1 (0 startに注意)
   
      call VecGetArrayF90(v,p_v,ierr)
      p_v(1) = (2d0, -4d0)
      call VecRestoreArrayF90(v,p_v,ierr)
   
      call VecView(v, PETSC_VIEWER_STDOUT_SELF, ierr)
   
      call VecDestroy(v,ierr)
      call PetscFinalize(ierr)
   
  end program petsc_test