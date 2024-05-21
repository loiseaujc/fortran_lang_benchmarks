module rhofunc
implicit none
public
integer, parameter :: dp=kind(0.d0)
contains
    pure real(dp) function rho(x,y)
        real(dp), intent(in) :: x,y
        if (x > 0.6_dp .and. x < 0.8_dp .and. y > 0.6_dp .and. y<0.8_dp) then
            rho = 1.0_dp
        else if (x> 0.2_dp .and. x<0.4_dp .and. y>0.2_dp .and. y<0.4_dp) then
            rho = -1.0_dp
        else
            rho = 0.0_dp
        end if
    end function

end module

program poisson
    use rhofunc, only: rho
    implicit none
    integer, parameter :: dp = kind(0.d0)

    !----- Physical parameters -----
    integer , parameter :: m = 300
    !! Number of grid points in each direction. 
    real(dp), parameter :: dx = 0.01_dp !1.0_dp / (m-1)
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

    ! Tic.
    call cpu_time(start_time)

    ! Initialize right-hand side.
    do j = 1, m
        do i = 1, m
            rhs(i,j) = rho(i*dx,j*dx)
        end do
    end do

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
    rhs = dx**2/4/epsilon0 * rhs

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

    write(*, *) "Solution converged to the desired tolerance after", iteration, "iterations."
    
    ! Toc.
    call cpu_time(end_time)

    write(*, *) "Computation required", end_time-start_time, "seconds (wall clock time)."

end program
