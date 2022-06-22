program main
#include <petsc/finclude/petscksp.h>
	use petscksp
	use user_module, only: User, SparseForm
	use solver_module, only: UserInitializeLinearSolver, UserDoLinearSolver, UserFinalizeLinearSolver
  implicit none

!    User-defined context that contains all the data structures used
!    in the linear solution process.

!   Vec    x,b      /* solution vector, right hand side vector and work vector */
!   Mat    A        /* sparse matrix */
!   KSP   ksp     /* linear solver context */
!   int    m,n      /* grid dimensions */
!
!   Since we cannot store Scalars and integers in the same context,
!   we store the integers/pointers in the user-defined context, and
!   the scalar values are carried in the common block.
!   The scalar values in this simplistic example could easily
!   be recalculated in each routine, where they are needed.
!
!   Scalar hx2,hy2  /* 1/(m+1)*(m+1) and 1/(n+1)*(n+1) */

!  Note: Any user-defined Fortran routines MUST be declared as external.

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                   Variable declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  PetscScalar  hx,hy,x,y
  type(User) userctx
  PetscErrorCode ierr
  PetscInt m, n, t, tmax, i, j
  PetscReal  enorm
  PetscScalar cnorm
  PetscScalar, ALLOCATABLE :: userx(:)
  PetscScalar, ALLOCATABLE :: userb(:)
  PetscScalar, ALLOCATABLE :: solution(:)
  type(SparseForm), allocatable :: form(:)

  PetscReal hx2, hy2
  common /param/ hx2, hy2

  tmax = 2
  N = 5

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  if (ierr /= 0) then
    print*,'Unable to initialize PETSc'
    stop
  endif

!  The next two lines are for testing only; these allow the user to
!  decide the grid size at runtime.

!  Create the empty sparse matrix and linear solver data structures

  call UserInitializeLinearSolver(N, userctx, ierr)
	CHKERRA(ierr)

!  Allocate arrays to hold the solution to the linear system.  This
!  approach is not normally done in PETSc programs, but in this case,
!  since we are calling these routines from a non-PETSc program, we
!  would like to reuse the data structures from another code. So in
!  the context of a larger application these would be provided by
!  other (non-PETSc) parts of the application code.

  allocate (userx(N), userb(N), solution(N), form(9))

!  Allocate an array to hold the coefficients in the elliptic operator


!  right-hand-side b[] and the solution with a known problem for testing.

  do i = 1, N
  !  userb(i) = 1.0d0
  !  solution(i) = 1.0d0
  !!  form(i)%ColumnId = i-1
  !  form(i)%RowId = i-1
  !  form(i)%Value = 1.0d0
  end do

  form(1)%ColumnId = 0; form(1)%RowId = 0; form(1)%Value = 1.0d0
  form(2)%ColumnId = 1; form(2)%RowId = 0; form(2)%Value = 5.0d0
  form(3)%ColumnId = 1; form(3)%RowId = 1; form(3)%Value = 2.0d0
  form(4)%ColumnId = 2; form(4)%RowId = 1; form(4)%Value = 8.0d0
  form(5)%ColumnId = 2; form(5)%RowId = 2; form(5)%Value = 3.0d0
  form(6)%ColumnId = 3; form(6)%RowId = 2; form(6)%Value = 9.0d0
  form(7)%ColumnId = 3; form(7)%RowId = 3; form(7)%Value = 4.0d0
  form(8)%ColumnId = 4; form(8)%RowId = 3; form(8)%Value = 10.0d0
  form(9)%ColumnId = 4; form(9)%RowId = 4; form(9)%Value = 5.0d0

  userb(1) = 11.0
  userb(2) = 28.0
  userb(3) = 45.0
  userb(4) = 66.0
  userb(5) = 25.0

  solution(1) = 1.0d0
  solution(2) = 2.0d0
  solution(3) = 3.0d0
  solution(4) = 4.0d0
  solution(5) = 5.0d0

!  Loop over a bunch of timesteps, setting up and solver the linear
!  system for each time-step.
!  Note that this loop is somewhat artificial. It is intended to
!  demonstrate how one may reuse the linear solvers in each time-step.

	!do t = 1, tmax
    call UserDoLinearSolver(form, userctx, userb, userx, ierr)
		CHKERRA(ierr)

!        Compute error: Note that this could (and usually should) all be done
!        using the PETSc vector operations. Here we demonstrate using more
!        standard programming practices to show how they may be mixed with
!        PETSc.
    cnorm = 0.0
    do i = 1, N
      cnorm = cnorm + PetscConj(solution(i)-userx(i))*(solution(i)-userx(i))
		enddo
    enorm =  PetscRealPart(cnorm*hx*hy)
    write(6,115) m,n,enorm
 115     format ('m = ',I2,' n = ',I2,' error norm = ',1PE11.4)
	!enddo

!  We are finished solving linear systems, so we clean up the
!  data structures.

  DEALLOCATE (userx,userb,solution)

  call UserFinalizeLinearSolver(userctx,ierr)
	CHKERRA(ierr)
  call PetscFinalize(ierr)
end program main
