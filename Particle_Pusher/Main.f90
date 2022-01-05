PROGRAM MAIN
  
  USE ISO_FORTRAN_ENV
  USE domain_tools
  USE GAUSS_SEIDEL
  USE write_netcdf
  USE VELOCITY_VERLET
  USE command_line ! script that enables the writing of variables on the command line
                   ! Written by Chris Brady and Heather Ratcliffe University of Warwick
  
  
  IMPLICIT NONE
  
  INTEGER :: n_x, n_y, gs_iter, ver_iter
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: E_x, E_y, phi, rho ! Gauss-Seidel Variables
  REAL(KIND=REAL64), DIMENSION(0:1000) :: x, y, v_x, v_y, a_x, a_y ! Verlet variables
  REAL(KIND=REAL64) :: d_t, d_x, d_y ! grid and time spacing
  CHARACTER(LEN=100) :: problem
  CHARACTER(LEN=*), PARAMETER :: filename = "particle_pusher_data.nc" ! filename for netcdf
  LOGICAL :: n_x_pass, n_y_pass, problem_pass
  REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: grid_X, grid_Y


  

  CALL parse_args

  ! conditonals for the values set using the command line
  n_x_pass = get_arg("n_x", n_x)
  n_y_pass = get_arg("n_y", n_y)
  problem_pass = get_arg("problem", problem)


  ! conditions for ensuring the passig and input of command line variables

  IF (n_x_pass) THEN
    PRINT*, "n_x=",n_x
  ELSE
    PRINT*, "PLEASE INPUT n_x"
    STOP "Stopped"
  END IF
  

  IF (n_y_pass) THEN
    PRINT*, "n_y=",n_y
  ELSE
    PRINT*, "PLEASE INPUT n_y"
    STOP "Stopped"
  END IF

  IF (problem_pass) THEN
    PRINT*, "problem=", problem
  ELSE
    PRINT*, "PLEASE INPUT PROBLEM AS 'null', 'single' or 'double'"
    STOP "Stopped"
  END IF


  ! allocate without ghost cells
  ALLOCATE(E_x(1:n_x, 1:n_y))
  ALLOCATE(E_y(1:n_x, 1:n_y))
  
  ! output charge density (rho) and potential (phi)

  CALL create_axis(grid_X, nels=n_x, axis_range=axis_range, nghosts=nghosts, delta=d_x)
  CALL create_axis(grid_Y, nels=n_y, axis_range=axis_range, nghosts=nghosts, delta=d_y)
  CALL get_charge_density(n_x, n_y, problem, rho)
  CALL get_potential(n_x, n_y, phi, problem, gs_iter)
  
  ! use functions to output E
  E_x = get_electric_field_x(n_x, n_y, phi)
  E_y = get_electric_field_y(n_x, n_y, phi)

  ! outputs verlet algorithm 
  CALL verlet_algorithm(n_x, n_y, E_x, E_y, problem, d_t, ver_iter, x, y, v_x, v_y, a_x, a_y)
  
  ! writes netcdf file 
  CALL writer_prototype(filename, rho, phi, E_x, E_y, x, y, v_x, v_y, a_x, a_y, n_x, n_y, problem, ver_iter, d_x, d_y, d_t)
  
  
  
  
END PROGRAM MAIN


