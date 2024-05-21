program poisson
    !! This program solves the Poisson equation on the unit square with
    !! homogeneous Dirichlet boundary conditions. The Laplace operator
    !! is discretized with the standard second-order accurate central
    !! finite difference scheme. The resulting discrete problem is solved
    !! using the Jacobi iterative solver. The numerical implementation
    !! relies on relatively standard Fortran constructs and is very similar
    !! to what would be asked from students in a numerical analysis class.
    implicit none
    !----- Sets the double precision kind -----
    integer, parameter :: precision = 15, range = 307
    integer, parameter :: dp = selected_real_kind(precision, range)
 
    !----- Physical parameters -----
    integer , parameter :: m = 300
    !! Number of grid points in each direction. 
    real(dp), parameter :: dx = 1.0_dp / (m-1)
    !! Uniform grid size in each direction.
    real(dp) :: rhs(m, m)
    !! Right-hand side of the Poisson equation.
    real(dp) :: phi(m, m)
    !! Solution of the Poisson equation.

    !----- Jacobi solver -----
    real(dp), parameter :: tolerance = 1e-6_dp
    !! Tolerance for the iterative Jacobi solver.
    real(dp) :: phi_prime(m, m)
    !! Working array.
    real(dp) :: residual
    !! Residual is defined as the infinity-norm of the update, i.e. max | phi - phi' |
    integer  :: iteration

    !----- Miscellaneous -----
    real(dp) :: start_time, end_time
    integer             :: i, j
    real(dp), parameter :: epsilon0 = 8.85E-12_dp

    ! Initialize right-hand side.
    do j = 1, m
        do i = 1, m
            rhs(i,j) = rho(i*dx,j*dx)
        end do
    end do

    ! Tic.
    call cpu_time(start_time)

    !---------------------------------
    !-----     JACOBI SOLVER     -----
    !---------------------------------

    ! Iteration counter.
    iteration = 0
    ! Working arrays.
    phi = 0.0_dp ; phi_prime = 0.0_dp
    ! Residual (ensuring at least one iteration is performed).
    residual = 1.0_dp
    ! Rescale rhs.
    rhs = dx**2/(4*epsilon0) * rhs

    ! Iterative solver.
    do while (residual > tolerance )
        ! Update iteration counter.
        iteration = iteration + 1
        ! Reset residual.
        residual = 0.0_dp
        ! Jacobi iteration.
        do j = 2, m-1
            do i = 2, m-1
                ! Jacobi update.
                phi_prime(i,j) = (phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1))/4 + rhs(i, j)
                ! On-the-fly residual computation.
                residual = max(residual, abs(phi_prime(i,j) - phi(i,j)))
            end do
        end do
        ! Updated solution.
        phi = phi_prime 
    end do
   
    ! Toc.
    call cpu_time(end_time)

    write(*, *) "Number of iterations  :", iteration
    write(*, *) "Wall clock time (sec) :", end_time-start_time

contains
     pure real(dp) function rho(x, y)
        implicit none
        real(dp), intent(in) :: x, y
        if (all([x, y] > 0.6_dp .and. [x, y] < 0.8_dp)) then
            rho = 1.0_dp
        else
            rho = merge(-1.0_dp, 0.0_dp, all([x, y]>0.2_dp .and. [x, y] < 0.4_dp))
        endif
   end function
end program
