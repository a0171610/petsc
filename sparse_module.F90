module sparse_module
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
  PetscErrorCode   ierr
contains

  subroutine solve_linear_equation(A, x, b, n)
    use petscksp
    integer, intent(in) :: n
    real(8), intent(in) :: A(n, n), b(n)
    real(8), intent(out) :: x(n)
    Mat :: At
    Vec :: xt, bt

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call transform_to_petsc_matrix(A, At, n)
    call transform_to_petsc_vector(b, bt, n)
    call transform_to_petsc_vector(x, xt, n)
    call solve(At, xt, bt)
  end subroutine solve_linear_equation

  subroutine transform_to_petsc_matrix(A, At, n)
    use petscksp
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n, n)
    Mat :: At
    integer :: i, idx(n)
  
  
    do i = 1, n
      idx(i) = i - 1
    enddo
  
    call MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, n, PETSC_NULL_INTEGER, At, ierr)
    do i = 1, n
      call MatSetValues(At, 1, i - 1, n, idx, A(i, :), INSERT_VALUES,ierr)
    enddo
  
    call MatAssemblyBegin(At, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(At, MAT_FINAL_ASSEMBLY, ierr)
  
  end subroutine transform_to_petsc_matrix

  subroutine transform_to_petsc_vector(x, xt, n)
    use petscksp
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    integer :: idx(n), i
    PetscScalar :: v(n)
    Vec :: xt
  
    do i = 1, n
      idx(i) = i - 1
      v(i) = x(i)
    enddo
  
    call VecCreateSeq(PETSC_COMM_SELF, n, xt, ierr)
    call VecSetValues(xt, n, idx, v, INSERT_VALUES,ierr)
  
  end subroutine transform_to_petsc_vector

  subroutine solve(At, xt, bt)
    use petscksp
    implicit none
    Mat :: At
    Vec :: xt, bt
    KSP :: ksp !連立方程式ソルバ構造体
    PC  :: pc  !前処理構造体
    KSPConvergedReason :: res
    integer ::iter

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KSP prep. (Krylov SubsPace method)
    call KSPCreate(PETSC_COMM_SELF, ksp, ierr)
    call KSPSetType(ksp, KSPBICG, ierr) !何も指定しないとGMRESが入る
    !KSPTYPEはKSPのあとに, CG, BICG, MINRES, GMRES, BCGS(=biCgStab)
 
    call KSPSetOperators(ksp, At, At, ierr) !前のaは解く行列,後ろは前処理用(同じでよい)
    call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr) !指定しないとゼロ初期値で解く
 
    call KSPGetPC(ksp, pc, ierr) 
    call PCSetType(pc, PCJACOBI, ierr) 
    !PCTYPEはPCのあとにJACOBI, ILU, ICC(不完全cholesky), LU, CHOLESKY(直接法)
 
    call KSPSetTolerances(ksp,1d-9,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
    call KSPSetFromOptions(ksp, ierr)
    call KSPView(ksp, PETSC_VIEWER_STDOUT_SELF, ierr) !オプションのレビュー
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Solve
    call KSPSolve(ksp, bt, xt, ierr)
 
    call KSPGetIterationNumber(ksp, iter, ierr)
    call KSPGetConvergedReason(ksp, res, ierr) !正なら収束, 負なら発散
    ! 1,2:rtol, 3,9:atol, -3:its, -4:dtol, -5,6:breakdown, -9:nan/inf
    ! -7:non-sym, -8,10:indefinite, -11:pcfail
 
    call MatView(At, PETSC_VIEWER_STDOUT_SELF, ierr)
    call VecView(xt, PETSC_VIEWER_STDOUT_SELF, ierr)
    call VecView(bt, PETSC_VIEWER_STDOUT_SELF, ierr)
  end subroutine solve
end module sparse_module