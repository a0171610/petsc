module user_module
#include <petsc/finclude/petscksp.h>
  use petscksp
  type User
    Vec x
    Vec b
    Mat A
    KSP ksp
    PetscInt N
  end type User

  type SparseForm
    integer :: RowId
    integer ::  ColumnId
    real(8) ::  Value
  end type SparseForm
end module user_module