program main
  use sparse_module, only: solve_linear_equation
  implicit none
  integer, parameter :: n = 3
  real(8) :: A(n, n), x(n), b(n)
  integer :: i

  A(1, 1) = 8.0d0
  A(1, 2) = 1.0d0
  A(1, 3) = 6.0d0
  A(2, 1) = 3.0d0
  A(2, 2) = 5.0d0
  A(2, 3) = 7.0d0
  A(3, 1) = 4.0d0
  A(3, 2) = 9.0d0
  A(3, 3) = 2.0d0

  b(:) = 1.0d0
  x(:) = 0.0d0

  call solve_linear_equation(A, x, b, n)

  do i = 1, 3
    write(*,*) x(i)
  end do
  
end program main