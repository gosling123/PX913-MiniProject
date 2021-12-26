MODULE GAUSS_SEIDEL !Module to handle claculation of phi, E etc.

    USE ISO_FORTRAN_ENV
    USE domain_tools
    IMPLICIT NONE
    SAVE
  
    REAL(KIND=REAL64), DIMENSION(2) :: axis_range = (/-1.0_REAL64, 1.0_REAL64/)
    INTEGER, PARAMETER :: nghosts = 1
    
    CONTAINS

    !Subroutine to genearte grid structure using create axis supplied by
    ! Chris Brady and Heather Ratcliffe (University of Warwick)
  
    SUBROUTINE get_grid(n_x, n_y, X, Y, delta_x, delta_y)
  
      INTEGER, INTENT(IN) :: n_x, n_y !Number of grid elemnts on axis
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y
      REAL(KIND=REAL64), INTENT(OUT) :: delta_x, delta_y !Grid spacing on each axis
  

      !Creates X-axis
      CALL create_axis(X, nels=n_x, axis_range=axis_range, nghosts=nghosts, delta=delta_x)
  
      !Creates Y-axis
      CALL create_axis(Y, nels=n_y, axis_range=axis_range, nghosts=nghosts, delta=delta_y)
  
  
    END SUBROUTINE get_grid
  
  
  
  
    !Subroutine to handle the geneartion of selcted charge densities
    SUBROUTINE get_charge_density(n_x, n_y, problem, rho)
  
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y
      INTEGER :: i, j
      INTEGER, INTENT(IN) :: n_x, n_y
      CHARACTER(LEN=*), INTENT(IN) :: problem !defines the set up of rho
      REAL(KIND=REAL64) :: exponent_0, exponent_1, exponent_2
      REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: rho
      REAL(KIND=REAL64) :: delta_x, delta_y 
        
      CALL get_grid(n_x, n_y, X, Y, delta_x, delta_y)
  
  
      ALLOCATE(rho(0:n_x+1, 0:n_y+1)) ! allocate with ghost cells indexed as 0 and n_x+1 or n_y+1 depending on axis
  
      ! null -- charge density is equal to zero in all space
      IF (problem .EQ. "null") THEN
  
        rho = 0.0_REAL64


      ! single -- chrage density has spaitial variation defined by a single Gaussian function
      ELSE IF (problem .EQ. "single") THEN
  
        rho = 0.0_REAL64
          
          
        DO j=1, n_y
          DO i=1, n_x
  
            exponent_0 = ((X(i)/0.1)**2 + (Y(j)/0.1)**2)
  
            rho(i,j) = EXP(-1.0_REAL64*exponent_0)
  
          END DO 
        END DO 
  
  
  
      !double -- chrage density has spaitial variation defined by a double Gaussian function  
      ELSE IF (problem .EQ. "double") THEN
  
        rho = 0.0_REAL64
  
        DO j=1, n_y
          DO i=1, n_x
  
            exponent_1 = -1.0_REAL64*(((X(i)+0.25)/0.1)**2 + ((Y(j)+0.25)/0.1)**2) 
            exponent_2 = -1.0_REAL64*(((X(i)-0.75)/0.2)**2 + ((Y(j)-0.75)/0.2)**2)
  
            rho(i,j) = EXP(exponent_1) + EXP(exponent_2)
          END DO 
        END DO
  
      ELSE
        ! Need warning to make sure user inputs correct options
        PRINT*, "WARNING: Please set problem as null, single or double"
  
      END IF
  
  
      
    END SUBROUTINE get_charge_density
  
  
    !subroutine to calculate electric potential from Poission's Equation using Gauss-Seidel iteration
    SUBROUTINE get_potential(n_x, n_y, phi, problem)
  
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y
      INTEGER :: i, j, k, N
      INTEGER, INTENT(IN) :: n_x, n_y
      CHARACTER(LEN=*), INTENT(IN) :: problem
      REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: rho
      REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: phi
      REAL(KIND=REAL64) :: delta_x, delta_y, x_part, y_part, denominator
      REAL(KIND=REAL64) :: e_tot, d_rms, x_FD, y_FD, x_error, y_error, d
  
      REAL(KIND=REAL64), PARAMETER :: tol = 0.00001_REAL64
      INTEGER, PARAMETER :: iterations = 1000000
  
      CALL get_charge_density(n_x, n_y, problem, rho) !required to claculate phi from Poissions Equations
  
      CALL get_grid(n_x, n_y, X, Y, delta_x, delta_y)
  
      ALLOCATE(phi(0:n_x+1, 0:n_y+1)) ! allocate with ghost cells indexed as 0 and n_x+1 or n_y+1 depending on axis
  
  
      phi = 0.0_REAL64

      !Denominator of Gauss-Seidel iteration
      denominator = 2/(delta_x)**2 + 2/(delta_y)**2  
       
      !e_tot - error defined as the sum of the differnces between the finite difference form of the laplacian of phi
      ! and the value of rho at each point on the grid.
      e_tot = 1.0_REAL64

      !d_rms - root mean square of the calculated error in phi
      d_rms = 1.0_REAL64
  
      N = n_x*n_y 
 
      !Start of Gauss-Seidel process
      DO k=1, iterations
  
        IF (ABS(e_tot/d_rms) <= tol) EXIT !condition for breaking loop and outputting a converged value of phi
  
        e_tot = 0.0_REAL64
        d_rms = 0.0_REAL64
    
        DO j=1, n_y
          DO i=1, n_x
              
            x_part = (phi(i+1,j) + phi(i-1,j))/(delta_x)**2
            y_part = (phi(i,j+1) + phi(i,j-1))/(delta_y)**2
  
            phi(i,j) = -1.0_REAL64*(rho(i,j)-x_part-y_part)/denominator !Gauss-Seidel Iteartion of phi
  
          END DO 
        END DO
  
  
        DO j=1, n_y
          DO i=1, n_x
  
            x_FD = (phi(i+1,j) - 2*phi(i,j) + phi(i-1,j))/(delta_x)**2 !second order x deriavtive             
            y_FD = (phi(i,j+1) - 2*phi(i,j) + phi(i,j-1))/(delta_y)**2 !second order y deriavtive
  
            x_error =  x_FD
            y_error =  y_FD
  
            e_tot = e_tot + ABS(x_error + y_error - rho(i,j))  ! summing to find total error
  
            d = SQRT(1.0_REAL64/REAL(N, REAL64))*SQRT(x_error**2 + y_error**2)
  
            d_rms = d_rms + d ! summing to find total root mean square
  
          END DO 
        END DO
  
        PRINT*,'Iteration:', k
        print*, 'e_tot', e_tot,'d_rms', d_rms, 'Condition', e_tot/d_rms
  
      
  
        ! Warning if phi has failed to converge
        IF(k >= iterations) THEN
          PRINT*, "WARNING, NOT CONVERGED"
        END IF
  
  
        
      END DO 
  
  
    END SUBROUTINE get_potential
  
    !function to caluate the electric field in x from phi taken as an argument from calling phi in
    !the main program 
    FUNCTION get_electric_field_x(n_x, n_y, phi)
       
      INTEGER, INTENT(IN) :: n_x, n_y
      REAL(KIND=REAL64), DIMENSION(0:n_x+1,0:n_y+1) :: phi
      REAL(KIND=REAL64), DIMENSION(1:n_x, 1:n_y):: E_x, get_electric_field_x !Alloacte in main program
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y 
      REAL(KIND=REAL64) :: delta_x, delta_y
      INTEGER :: i, j
  
  
      CALL get_grid(n_x, n_y, X, Y, delta_x, delta_y)
  
      DO j=1, n_y
        DO i=1, n_x
  
          E_x(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*delta_x) !Electric field in x
  
        END DO 
      END DO
  
      get_electric_field_x = E_x
  
    END FUNCTION get_electric_field_x

    !function to caluate the electric field in x from phi taken as an argument from calling phi in
    !the main program 
    FUNCTION get_electric_field_y(n_x, n_y, phi)
       
      INTEGER, INTENT(IN) :: n_x, n_y
      REAL(KIND=REAL64), DIMENSION(0:n_x+1,0:n_y+1) :: phi
      REAL(KIND=REAL64), DIMENSION(1:n_x, 1:n_y) :: E_y, get_electric_field_y !Alloacte in main program
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y 
      REAL(KIND=REAL64) :: delta_x, delta_y
      INTEGER :: i, j
  
      CALL get_grid(n_x, n_y, X, Y, delta_x, delta_y) 
  
      DO i=1, n_x
        DO j=1, n_y
  
          E_y(i,j) = (phi(i,j+1) - phi(i,j-1))/(2*delta_y) !Electric field in y
  
        END DO 
      END DO
  
      get_electric_field_y = E_y
  
    END FUNCTION get_electric_field_y
  
    
  END MODULE GAUSS_SEIDEL
  
    
  
  
  
