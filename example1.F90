program petsc_mat
#include <petsc/finclude/petscksp.h>
      use petscksp
implicit none 
 
    integer, parameter :: n = 3
    integer :: ierr, iter
    PetscScalar :: v(n), zero 
    Mat :: a
    Vec :: b, x
 
    KSP :: ksp !連立方程式ソルバ構造体
    PC  :: pc  !前処理構造体
    KSPConvergedReason :: res
 
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Mat prep.
    call MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,n,PETSC_NULL_INTEGER,a,ierr)
    v(:) = (/8d0, 1d0, 6d0/)
    call MatSetValues(a, 1, 0, n, (/0,1,2/), v, INSERT_VALUES,ierr)
    v(:) = (/3d0, 5d0, 7d0/)
    call MatSetValues(a, 1, 1, n, (/0,1,2/), v, INSERT_VALUES,ierr)
    v(:) = (/4d0, 9d0, 2d0/)
    call MatSetValues(a, 1, 2, n, (/0,1,2/), v, INSERT_VALUES,ierr)
    call MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(a,MAT_FINAL_ASSEMBLY,ierr)
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Vec prep.
    call VecCreateSeq(PETSC_COMM_SELF, n, b,ierr)
    v(:) = (/1d0, 1d0, 1d0/)
    call VecSetValues(b, n, (/0,1,2/), v, INSERT_VALUES,ierr)
 
    call VecDuplicate(b,x,ierr)
    zero = 0d0 !1d0/15d0 !で初期値ありにするとイタレーションが0になる
    call VecSet(x,zero,ierr) !あらかじめ初期値を指定できる
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KSP prep. (Krylov SubsPace method)
    call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
    call KSPSetType(ksp, KSPBICG, ierr) !何も指定しないとGMRESが入る
    !KSPTYPEはKSPのあとに, CG, BICG, MINRES, GMRES, BCGS(=biCgStab)
 
    call KSPSetOperators(ksp,a,a,ierr) !前のaは解く行列,後ろは前処理用(同じでよい)
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr) !指定しないとゼロ初期値で解く
 
    call KSPGetPC(ksp,pc,ierr) 
    call PCSetType(pc,PCJACOBI,ierr) 
    !PCTYPEはPCのあとにJACOBI, ILU, ICC(不完全cholesky), LU, CHOLESKY(直接法)
 
    call KSPSetTolerances(ksp,1d-9,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
    call KSPSetFromOptions(ksp,ierr)
    call KSPView(ksp,PETSC_VIEWER_STDOUT_SELF,ierr) !オプションのレビュー
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Solve
    call KSPSolve(ksp,b,x,ierr)
 
    call KSPGetIterationNumber(ksp,iter,ierr)
    call KSPGetConvergedReason(ksp, res, ierr) !正なら収束, 負なら発散
    ! 1,2:rtol, 3,9:atol, -3:its, -4:dtol, -5,6:breakdown, -9:nan/inf
    ! -7:non-sym, -8,10:indefinite, -11:pcfail
 
    call MatView(a, PETSC_VIEWER_STDOUT_SELF, ierr)
    call VecView(x, PETSC_VIEWER_STDOUT_SELF, ierr)
    call VecView(b, PETSC_VIEWER_STDOUT_SELF, ierr)
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Finalize
    call KSPDestroy(ksp,ierr)
    call MatDestroy(a,ierr)
    call VecDestroy(b,ierr)
    call PetscFinalize(ierr)
 
end program petsc_mat