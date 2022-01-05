MODULE VELOCITY_VERLET !Module to handle claculation of phi, E etc.

  USE ISO_FORTRAN_ENV
  USE GAUSS_SEIDEL
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: m = 1, q = -1
  
  CONTAINS

  
  !subroutine to define initial velovity and posistions
  SUBROUTINE initial_condition(x_0, v_0, problem)

    CHARACTER(LEN=*), INTENT(IN) :: problem !defines the set up of problem
    REAL(KIND=REAL64), DIMENSION(2), INTENT(OUT) :: x_0, v_0 !initial position and velocity

    IF (problem .EQ. "null") THEN ! for rho = 0 

      x_0 = (/0.0_REAL64, 0.0_REAL64/)
      v_0 = (/0.1_REAL64, 0.1_REAL64/)

    ELSE IF (problem .EQ. "single") THEN ! for rho with single gaussian

      x_0 = (/0.1_REAL64, 0.0_REAL64/)
      v_0 = (/0.0_REAL64, 0.0_REAL64/)

  
    ELSE IF (problem .EQ. "double") THEN ! for rho with double gaussian

      x_0 = (/0.0_REAL64, 0.5_REAL64/)
      v_0 = (/0.0_REAL64, 0.0_REAL64/)

    END IF


  END SUBROUTINE initial_condition


  !subroutine to perform the velocity verlet algorithm
  SUBROUTINE verlet_algorithm(n_x, n_y, E_x, E_y, problem, d_t, ver_iter, x, y, v_x, v_y, a_x, a_y)
     
      INTEGER, PARAMETER :: max_iter = 1000 !set number of iterations
      REAL(KIND=REAL64), INTENT(OUT):: d_t !time step
      
      INTEGER, INTENT(IN) :: n_x, n_y
      CHARACTER(LEN=*), INTENT(IN) :: problem !defines the set up of problem
      REAL(KIND=REAL64), DIMENSION(1:n_x, 1:n_y), INTENT(IN) :: E_x, E_y
      REAL(KIND=REAL64), DIMENSION(0:1000), INTENT(OUT) :: x, y !posistions
      REAL(KIND=REAL64), DIMENSION(0:1000), INTENT(OUT) :: v_x, v_y  !velocity
      REAL(KIND=REAL64), DIMENSION(0:1000), INTENT(OUT) :: a_x, a_y !acceleration
      INTEGER :: cell_x, cell_y, k 
      INTEGER,INTENT(OUT) :: ver_iter !number of iterations completed
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: grid_X, grid_Y !grid 
      REAL(KIND=REAL64), DIMENSION(2) :: x_0, v_0 !initial conditions
      REAL(KIND=REAL64) :: d_x, d_y !grid spacing

      !initialsise all elements as zero
      x = 0.0_REAL64
      v_x = 0.0_REAL64
      a_x = 0.0_REAL64

      y = 0.0_REAL64
      v_y = 0.0_REAL64
      a_y = 0.0_REAL64

      !set time step
      d_t = 0.01_REAL64

      !get initial conditions
      CALL initial_condition(x_0, v_0, problem)
      
      !Creates grid_X-axis
      CALL create_axis(grid_X, nels=n_x, axis_range=axis_range, nghosts=nghosts, delta=d_x)

      !Creates grid_Y-axis
      CALL create_axis(grid_Y, nels=n_y, axis_range=axis_range, nghosts=nghosts, delta=d_y)


      !find indicies on the 2d grid
      cell_x = FLOOR((x_0(1)+1.0_REAL64)/d_x + 1)
      cell_y = FLOOR((x_0(2)+1.0_REAL64)/d_y + 1)
      
      !initial acceleration
      a_x(0) = (q/m)*E_x(cell_x, cell_y)
      a_y(0) = (q/m)*E_y(cell_x, cell_y)

      !initial posisition
      x(0) = x_0(1)
      y(0) = x_0(2)

      !initial velocity
      v_x(0) = v_0(1)
      v_y(0) = v_0(2)

      ver_iter = 0

      DO k = 1, max_iter
        ver_iter = ver_iter + 1 !closk of number of iterations performed

        !verlet update on posistions
        x(k) = x(k-1) + v_x(k-1)*d_t + 0.50_REAL64*a_x(k-1)*d_t**2.0_REAL64
        y(k) = y(k-1) + v_x(k-1)*d_t + 0.50_REAL64*a_y(k-1)*d_t**2.0_REAL64

        !find grid indicies
        cell_x = FLOOR((x(k)+1.0_REAL64)/d_x + 1)
        cell_y = FLOOR((y(k)+1.0_REAL64)/d_y + 1)

        ! exit loop if out of bounds
        IF (cell_x > n_x) EXIT
        IF (cell_y > n_y) EXIT

        ! PRINT*, cell_x_1, cell_x_2

        !verlet update on accelerations
        a_x(k) = (q/m)*E_x(cell_x, cell_y)
        a_y(k) = (q/m)*E_y(cell_x, cell_y)

        !verlet update on velocities
        v_x(k) = v_x(k-1) + 0.5_REAL64*d_t*(a_x(k) + a_x(k-1))
        v_y(k) = v_y(k-1) + 0.5_REAL64*d_t*(a_y(k) + a_y(k-1))


      END DO

      print*, ver_iter




  
  
    END SUBROUTINE verlet_algorithm







  
END MODULE VELOCITY_VERLET

  



