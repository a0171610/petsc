module solver_module
#include <petsc/finclude/petscksp.h>
  use petscksp
  use user_module, only: User
  implicit none

  PetscReal hx2, hy2

contains

  subroutine UserInitializeLinearSolver(N, userctx, ierr)
#include <petsc/finclude/petscksp.h>
    use user_module, only: User
    implicit none
  
    PetscInt N
    PetscErrorCode ierr
    type(User) userctx
    Mat A
    Vec b, x
    KSP ksp
  
    ! Local variable declararions
    PetscInt Ntot, five, one
  
    ! Here we assume use of a grid of size m x n, with all points on the
    ! interior of the domain, i.e., we do not include the points corresponding
    ! to homogeneous Dirichlet boundary conditions.  We assume that the domain
    ! is [0,1]x[0,1].
  
    Ntot = N
  
    five = 5
    one = 1
  
    ! Create the sparse matrix. Preallocate 5 nonzeros per row.

    ! 行列Aは Ntot * Ntot で、各行の非零要素がfive (多分最大数でok)
    call MatCreateSeqAIJ(PETSC_COMM_SELF,Ntot,Ntot,five,PETSC_NULL_INTEGER,A,ierr)
    CHKERRQ(ierr)
    !
    ! Create vectors. Here we create vectors with no memory allocated.
    ! This way, we can use the data structures already in the program
    ! by using VecPlaceArray() subroutine at a later stage.
    !
    call VecCreateSeqWithArray(PETSC_COMM_SELF,one,Ntot,PETSC_NULL_SCALAR,b,ierr)
    CHKERRQ(ierr)
    call VecDuplicate(b,x,ierr);CHKERRQ(ierr)
  
     ! Create linear solver context. This will be used repeatedly for all
    ! the linear solves needed.
  
    call KSPCreate(PETSC_COMM_SELF,ksp,ierr);CHKERRQ(ierr)
  
    userctx%x = x
    userctx%b = b
    userctx%A = A
    userctx%ksp = ksp
    userctx%N = N

  end subroutine UserInitializeLinearSolver
  ! -----------------------------------------------------------------------
  
  !   Solves -div (rho grad psi) = F using finite differences.
  !   rho is a 2-dimensional array of size m by n, stored in Fortran
  !   style by columns. userb is a standard one-dimensional array.
  
  subroutine UserDoLinearSolver(form, userctx, userb, userx, ierr)
#include <petsc/finclude/petscksp.h>
    use user_module, only: User, SparseForm
    implicit none
  
    PetscErrorCode ierr
    type(User) userctx
    type(SparseForm), allocatable :: form(:)
    PetscScalar userb(*),userx(*)
    Mat A
    Vec b,x
    KSP ksp
  
    PC   pc
    PetscInt N, one
    PetscInt i, j, row, col
    PetscScalar  val
    integer :: sz

    one  = 1
    x    = userctx%x
    b    = userctx%b
    A    = userctx%A
    ksp  = userctx%ksp
    N    = userctx%N
    sz = size(form)
 
  !  This is not the most efficient way of generating the matrix,
  !  but let's not worry about it.  We should have separate code for
  !  the four corners, each edge and then the interior. Then we won't
  !  have the slow if-tests inside the loop.
  !
  !  Compute the operator
  !          -div rho grad
  !  on an m by n grid with zero Dirichlet boundary conditions. The rho
  !  is assumed to be given on the same grid as the finite difference
  !  stencil is applied.  For a staggered grid, one would have to change
  !  things slightly.
  
    do i = 1, sz
      row = form(i)%RowId
      col = form(i)%ColumnId
      val = form(i)%Value
      call MatSetValues(A, one, row, one, col, val, INSERT_VALUES, ierr)
      CHKERRQ(ierr)
    enddo

  !
  !     Assemble matrix
  !
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  !
  !     Set operators. Here the matrix that defines the linear system
  !     also serves as the preconditioning matrix. Since all the matrices
  !     will have the same nonzero pattern here, we indicate this so the
  !     linear solvers can take advantage of this.
  !
    call KSPSetOperators(ksp,A,A,ierr);CHKERRQ(ierr)
  !
  !     Set linear solver defaults for this problem (optional).
  !     - Here we set it to use direct LU factorization for the solution
  !
    call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)
    call PCSetType(pc,PCLU,ierr);CHKERRQ(ierr)
  
  !
  !     Set runtime options, e.g.,
  !        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  !     These options will override those specified above as long as
  !     KSPSetFromOptions() is called _after_ any other customization
  !     routines.
  !
  !     Run the program with the option -help to see all the possible
  !     linear solver options.
  !
    call KSPSetFromOptions(ksp,ierr);CHKERRQ(ierr)
  
  !
  !     This allows the PETSc linear solvers to compute the solution
  !     directly in the user's array rather than in the PETSc vector.
  !
  !     This is essentially a hack and not highly recommend unless you
  !     are quite comfortable with using PETSc. In general, users should
  !     write their entire application using PETSc vectors rather than
  !     arrays.
  !
    call VecPlaceArray(x,userx,ierr);CHKERRQ(ierr)
    call VecPlaceArray(b,userb,ierr);CHKERRQ(ierr)
  
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !                      Solve the linear system
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call KSPSolve(ksp,b,x,ierr);CHKERRQ(ierr)

    call MatView(A, PETSC_VIEWER_STDOUT_SELF, ierr)
    call VecView(x, PETSC_VIEWER_STDOUT_SELF, ierr)
    call VecView(b, PETSC_VIEWER_STDOUT_SELF, ierr)
  
    call VecResetArray(x,ierr);CHKERRQ(ierr)
    call VecResetArray(b,ierr);CHKERRQ(ierr)
  end subroutine UserDoLinearSolver
  
  ! ------------------------------------------------------------------------
  
  subroutine UserFinalizeLinearSolver(userctx,ierr)
#include <petsc/finclude/petscksp.h>
    use user_module, only: User
    implicit none
  
    PetscErrorCode ierr
    type(User) userctx
  
  !
  !     We are all done and don't need to solve any more linear systems, so
  !     we free the work space.  All PETSc objects should be destroyed when
  !     they are no longer needed.
  !
    call VecDestroy(userctx%x,ierr);CHKERRQ(ierr)
    call VecDestroy(userctx%b,ierr);CHKERRQ(ierr)
    call MatDestroy(userctx%A,ierr);CHKERRQ(ierr)
    call KSPDestroy(userctx%ksp,ierr);CHKERRQ(ierr)
  end subroutine UserFinalizeLinearSolver
  
end module solver_module