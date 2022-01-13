module shared_data

    use iso_fortran_env
    implicit none

    !Global Constants
    real(kind=real64), parameter :: q = -1.0d0 !Charge
    real(kind=real64), parameter :: m =  1.0d0 !Mass
    real(kind=real64), parameter :: eps =  1.0d-16 !eps



    !System Variables
    integer :: n_x, n_y !Dimensions of Grid
    real(kind=real64) :: d_x, d_y !Grid Spacing
    real(kind=real64), allocatable :: grid_X(:), grid_Y(:) !Grid
    real(kind=real64), allocatable :: rho(:, :) !ndimsCharge Density, shape (0:n_x+1,0:n_y+1)
    real(kind=real64), allocatable :: phi(:, :) !Scalar Potential,shape (0:n_x+1,0:n_y+1)
    real(kind=real64), allocatable :: E(:,:, :) !Electric Field, shape(n_x, n_y, 2) 



    !Particle Trajectory Varaibles
    integer :: max_steps !Max number of steps
    integer :: viter ! Time step trajectory ends (due to leaving system)
    real(kind=real64) :: d_t=0.01_real64 ! Time step size, its set to -1 for error checking in setup_verlet
    real(kind=real64), allocatable :: r(:, :) ! Position, shape (max_iter, 2)
    real(kind=real64), allocatable :: v(:, :) !Velocity, shape (max_iter, 2)
    real(kind=real64), allocatable :: a(:, :) !Acceleration, shape (max_iter, 2)


end module 