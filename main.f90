program main
    
    use shared_data
    use VELOCITY_VERLET
    use GAUSS_SEIDEL
    use io
    use command_line
    implicit NONE

    character(len = 20) :: problem
    logical :: n_x_pass, n_y_pass, problem_pass
    CALL parse_args

    ! conditonals for the values set using the command line
    n_x_pass = get_arg("n_x", n_x)
    n_y_pass = get_arg("n_y", n_y)
    problem_pass = get_arg("problem", problem)
  
  
    ! conditions for ensuring the passing and input of command line variables
    IF (n_x_pass .or. n_x .le. 0) THEN
      PRINT*, "n_x=",n_x
    ELSE
      PRINT*, "PLEASE INPUT n_x (POSITIVE INTEGER)"
      STOP ""
    END IF
    
  
    IF (n_y_pass .or. n_y .le. 0) THEN
      PRINT*, "n_y=",n_y
    ELSE
      PRINT*, "PLEASE INPUT n_y (POSITIVE INTEGER)"
      STOP ""
    END IF
  
    IF (problem_pass) THEN
      PRINT*, "problem=", problem
    ELSE
      PRINT*, "PLEASE INPUT PROBLEM AS 'null', 'single' or 'double'"
      STOP ""
    END IF
  
    print *, "Program Started"


    call setup_sys()
    print *, "System Set up"

    call get_charge_density(problem)
    print *, "Rho Calculated"

    call get_potential()
    print *, "Potential Calculated"

    call get_electric_field()
    print *, "Electric Field Calculated"

    call setup_verlet(problem, 1000, 0.01d0)

    call iterate_verlet(0, 1000)

    call writer_prototype("particle_pusher_data.nc", problem)

end program