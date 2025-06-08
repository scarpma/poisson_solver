program poisson
  use precision
  implicit none
  real(dp), dimension(:,:), allocatable :: U, Uold
  real(dp), dimension(:,:), allocatable :: rho
  real(dp) :: tolerance, a, err, maxerr, w
  integer :: i,j,k, N, error
  character(1) :: M, BC
  character(20) :: arg
  real(dp), parameter :: Pi=3.141593_dp
  real(dp), parameter :: e2=14.4_dp ! eV*Ang
  real(dp) :: L, dx, D
  integer :: istart, iend, j1, j2
  
  ! command line arguments help
  if (iargc()<5) then
     write(*,*) 'poisson Method N tolerance w BC'
     write(*,*) '  - Method: J|G (Jacobi or Gauss-Siedel)'
     write(*,*) '  - N: <int> (number of points in x and y)'
     write(*,*) '  - tolerance: <real> (for convergence)'
     write(*,*) '  - w: <real> (relaxation only for G)'
     write(*,*) '  - BC: D|N (Dirichlet or Neumann)'
     stop
  endif
  
  ! parsing of command line arguments
  call getarg(1,arg)
  read(arg,*) M

  call getarg(2,arg)
  read(arg,*) N

  call getarg(3,arg)
  read(arg,*) tolerance

  call getarg(4,arg)
  read(arg,*) w

  call getarg(5,arg)
  read(arg,*) BC

  allocate(U(N,N), stat=error)
  allocate(Uold(N,N), stat=error)
  if (error /= 0) then
     write(*,*) 'allocation error'
     stop
  endif

  ! grid spacing
  dx = 1.0_dp / N

  ! initial condition for first iteration
  U=0.0_dp

  ! dirichlet boundary conditions (box)
  U(1:N,1) = 0.0_dp ! bottom edge
  U(1:N,N) = 0.0_dp ! top edge
  U(1,1:N) = 0.0_dp ! left edge
  U(N,1:N) = 0.0_dp ! right edge
  ! dirichlet boundary condition on plates
  L = 0.3_dp ! plate length
  D = 0.5_dp ! plate separation
  istart = int(N/2 - L/(2*dx))     ! Plate x-start
  iend   = int(N/2 + L/(2*dx))     ! Plate x-end
  j1     = int(N/2 - D/(2*dx))     ! Bottom plate y-index
  j2     = int(N/2 + D/(2*dx))     ! Top plate y-index
  U(istart:iend, j1) = -1.0_dp
  U(istart:iend, j2) = +1.0_dp

  select case (M)
  case("J")
    write(*,*) "Jacobi iteration"
  case("G")
    write(*,*) "Gauss-Siedel iteration"
  end select
  
  ! -----------------------------
  ! Iterative solver main loop
  ! -----------------------------
  err = 2.0_dp * tolerance
  k = 1
  do while (err > tolerance)
    Uold = U          ! Store previous solution
    maxerr = 0.0_dp   ! Reset max error
 
    ! ---------------
    ! loop in space
    ! ---------------
    do j = 2, N-1
      do i = 2, N-1

        ! Skip capacitor plates
        if (i >= istart .and. i<=iend .and. (j==j1 .or. j==j2)) cycle  
        
        ! solution domain: perform iteration update
        select case (M)
        case("J")
          U(i,j) = (Uold(i-1,j) + Uold(i+1,j) + Uold(i,j-1) + Uold(i,j+1))/4.0_dp
        case("G")
          U(i,j) = (U(i-1,j) + Uold(i+1,j) + U(i,j-1) + Uold(i,j+1))/4.0_dp
        end select

        ! relaxation factor w
        U(i,j) = (1-w)*Uold(i,j) + w*U(i,j)
        
        ! check covergence
        if (abs(Uold(i,j)-U(i,j)) > maxerr ) then
           maxerr = abs(Uold(i,j) - U(i,j))
        end if

        ! optional Neumann o(a^2) along x==1 and x==n
        if (BC.eq."N") then
          U(1,j) = 4.0_dp/3.0_dp * U(2,j) - 1.0_dp/3.0_dp * U(3,j)
          U(N,j) = 4.0_dp/3.0_dp * U(N-1,j) - 1.0_dp/3.0_dp * U(N-2,j)
          U(1,1) = 4.0_dp/3.0_dp * U(2,1) - 1.0_dp/3.0_dp * U(3,1)
          U(N,N) = 4.0_dp/3.0_dp * U(N-1,N) - 1.0_dp/3.0_dp * U(N-2,N)
        endif
   
      end do

      ! optional Neumann o(a^2) along y==1 and y==n
      if (BC.eq."N") then
        U(1:N/2-11,1) = 4.0_dp/3.0_dp * U(1:N/2-11,2) - 1.0_dp/3.0_dp * U(1:N/2-11,3)
        U(N/2+11:N,1) = 4.0_dp/3.0_dp * U(N/2+11:N,2) - 1.0_dp/3.0_dp * U(N/2+11:N,3)
        U(1:N/2-11,N) = 4.0_dp/3.0_dp * U(1:N/2-11,N-1) - 1.0_dp/3.0_dp * U(1:N/2-11,N-2)
        U(N/2+11:N,N) = 4.0_dp/3.0_dp * U(N/2+11:N,N-1) - 1.0_dp/3.0_dp * U(N/2+11:N,N-2)
      endif
    end do
 
    ! Output iteration progress
    write(*,*) 'iter: ',k, maxerr
    err = maxerr
    k = k + 1
  end do   
 





  ! write to file in fortran order
  open(101, file='sol.dat')
  do j = 1, N
     do i = 1, N
        write(101, *) U(i,j)
     end do
     write(101,*)
  end do
  close(101)
 
end program poisson
