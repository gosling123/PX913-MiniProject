module shared_data

    use iso_fortran_env
    implicit none


    !System Variables
    integer :: n_x, n_y !Dimensions of Grid
    real(kind=real64) :: d_x, d_y !Grid Spacing
    real(kind=real64), allocatable :: grid_X(:), grid_Y(:) !Grid
    real(kind=real64), allocatable :: rho(:, :) !ndimsCharge Density, shape (0:n_x+1,0:n_y+1)
    real(kind=real64), allocatable :: phi(:, :) !Scalar Potential,shape (0:n_x+1,0:n_y+1)
    real(kind=real64), allocatable :: E(:,:, :) !Electric Field, shape(n_x, n_y, 2) 



    !Particle Trajectory Varaibles
    integer :: max_steps !Max number of steps
    integer :: tstep_end ! Time step trajectory ends (due to leaving system)
    real(kind=real64) :: dt=-1.0_real64 ! Time step size, its set to -1 for error checking in setup_verlet
    real(kind=real64), allocatable :: r(:, :) ! Position, shape (max_iter, 2)
    real(kind=real64), allocatable :: v(:, :) !Velocity, shape (max_iter, 2)
    real(kind=real64), allocatable :: a(:, :) !Acceleration, shape (max_iter, 2)


    





end module 