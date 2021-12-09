MODULE GAUSS_SEIDEL

    USE ISO_FORTRAN_ENV
    USE domain_tools
    IMPLICIT NONE
    SAVE
  
    REAL(KIND=REAL64), DIMENSION(2) :: axis_range = (/-1.0_REAL64, 1.0_REAL64/)
    INTEGER, PARAMETER :: nghosts = 1
    
    CONTAINS
  
    SUBROUTINE get_grid(n_x, n_y, X, Y, delta_x, delta_y)
  
      INTEGER, INTENT(IN) :: n_x, n_y
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y
      REAL(KIND=REAL64), INTENT(OUT) :: delta_x, delta_y
  
      CALL create_axis(X, nels=n_x, axis_range=axis_range, nghosts=nghosts, delta=delta_x)
  
      CALL create_axis(Y, nels=n_y, axis_range=axis_range, nghosts=nghosts, delta=delta_y)
  
  
  
  
    END SUBROUTINE get_grid
  
  
  
  
  
    SUBROUTINE get_charge_density(n_x, n_y, problem, rho)
  
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y
      INTEGER :: i, j
      INTEGER, INTENT(IN) :: n_x, n_y
      CHARACTER(LEN=*), INTENT(IN) :: problem
      REAL(KIND=REAL64) :: exponent_0, exponent_1, exponent_2
      REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: rho
      REAL(KIND=REAL64) :: delta_x, delta_y 
        
      CALL get_grid(n_x, n_y, X, Y, delta_x, delta_y)
  
  
      ALLOCATE(rho(0:n_x+1, 0:n_y+1))
  
  
      IF (problem .EQ. "null") THEN
  
        rho = 0.0_REAL64
  
      ELSE IF (problem .EQ. "single") THEN
  
        rho = 0.0_REAL64
          
          
        DO i=1, n_x
          DO j=1, n_y
  
            exponent_0 = ((X(i)/0.1)**2 + (Y(j)/0.1)**2)
  
            rho(i,j) = EXP(-1.0_REAL64*exponent_0)
  
          END DO 
        END DO 
  
  
  
  
      ELSE IF (problem .EQ. "double") THEN
  
        rho = 0.0_REAL64
  
        DO i=1, n_x
          DO j=1, n_y
  
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
  
  
  
    SUBROUTINE get_potential(n_x, n_y, phi, problem)
  
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y
      INTEGER :: i, j, k, iterations, N
      INTEGER, INTENT(IN) :: n_x, n_y
      CHARACTER(LEN=*), INTENT(IN) :: problem
      REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: rho
      REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: phi
      REAL(KIND=REAL64) :: delta_x, delta_y, x_part, y_part, denominator
      REAL(KIND=REAL64) :: e_tot, d_rms, x_FD, y_FD, x_error, y_error, d
  
      REAL(KIND=REAL64), PARAMETER :: tol = 0.00001_REAL64
  
      CALL get_charge_density(n_x, n_y, problem, rho)
  
      CALL get_grid(n_x, n_y, X, Y, delta_x, delta_y)
  
      ALLOCATE(phi(0:n_x+1, 0:n_y+1))
  
  
      phi = 0.0_REAL64
  
      denominator = 2/(delta_x)**2 + 2/(delta_y)**2
  
      iterations = 1000000
  
      ! x_sum = 0
      ! y_sum = 0
      ! rho_sum = 0
  
      ! error = 0
  
      e_tot = 1.0_REAL64
      d_rms = 1.0_REAL64
  
      N = n_x*n_y
  
      DO k=1, iterations
  
        IF (ABS(e_tot/d_rms) <= tol) EXIT
  
        e_tot = 0.0_REAL64
        d_rms = 0.0_REAL64
  
        print*, e_tot, d_rms
    
        DO i=1, n_x
          DO j=1, n_y
              
            x_part = (phi(i+1,j) + phi(i-1,j))/(delta_x)**2
            y_part = (phi(i,j+1) + phi(i,j-1))/(delta_y)**2
  
            phi(i,j) = -1.0_REAL64*(rho(i,j)-x_part-y_part)/denominator
  
          END DO 
        END DO
  
  
        DO i=1, n_x
          DO j=1, n_y
  
            x_FD = (phi(i+1,j) - 2*phi(i,j) + phi(i-1,j))/(delta_x)**2
              
            y_FD = (phi(i,j+1) - 2*phi(i,j) + phi(i,j-1))/(delta_y)**2
  
            x_error =  x_FD
            y_error =  y_FD
  
            e_tot = e_tot + ABS(x_error + y_error - rho(i,j))   
  
            d = SQRT(1.0_REAL64/REAL(N, REAL64))*SQRT(x_error**2 + y_error**2)
  
            d_rms = d_rms + d
  
          END DO 
        END DO
  
  
        print*, e_tot, d_rms
  
        PRINT*, k
  
        
        ! print*, e_tot, d_rms
  
        
  
      
  
  
        IF(k >= iterations) THEN
          PRINT*, "WARNING, NOT CONVERGED"
        END IF
  
  
        
      END DO 
  
  
    END SUBROUTINE get_potential
  
  
    FUNCTION get_electric_field_x(n_x, n_y, phi)
       
      INTEGER, INTENT(IN) :: n_x, n_y
      REAL(KIND=REAL64), DIMENSION(0:n_x+1,0:n_y+1) :: phi
      REAL(KIND=REAL64), DIMENSION(1:n_x, 1:n_y):: E_x, get_electric_field_x
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y 
      REAL(KIND=REAL64) :: delta_x, delta_y
      INTEGER :: i, j
  
  
      CALL get_grid(n_x, n_y, X, Y, delta_x, delta_y)
  
      DO i=1, n_x
        DO j=1, n_y
  
          E_x(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*delta_x)
  
        END DO 
      END DO
  
      get_electric_field_x = E_x
  
    END FUNCTION get_electric_field_x
  
    FUNCTION get_electric_field_y(n_x, n_y, phi)
       
      INTEGER, INTENT(IN) :: n_x, n_y
      REAL(KIND=REAL64), DIMENSION(0:n_x+1,0:n_y+1) :: phi
      REAL(KIND=REAL64), DIMENSION(1:n_x, 1:n_y) :: E_y, get_electric_field_y
      REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE :: X, Y 
      REAL(KIND=REAL64) :: delta_x, delta_y
      INTEGER :: i, j
  
      CALL get_grid(n_x, n_y, X, Y, delta_x, delta_y)
  
      DO i=1, n_x
        DO j=1, n_y
  
          E_y(i,j) = (phi(i,j+1) - phi(i,j-1))/(2*delta_y)
  
        END DO 
      END DO
  
      get_electric_field_y = E_y
  
    END FUNCTION get_electric_field_y
  
    
  END MODULE GAUSS_SEIDEL
  
  


! used chris and heathers netcdf example for quick plot check
  MODULE write_netcdf
  
    USE ISO_FORTRAN_ENV
    USE netcdf
  
    IMPLICIT NONE
  
    CONTAINS
  
    SUBROUTINE writer_prototype(larr, filename, ierr)
  
      REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: larr
      INTEGER, PARAMETER :: ndims = 2
      ! We can use this parameter here, which makes it easier to 
      ! replicate this function for different dimensionalities
      INTEGER, DIMENSION(ndims) :: sizes, dim_ids
      CHARACTER(LEN=1), DIMENSION(ndims) :: dims=(/"x", "y" /)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      INTEGER :: ierr, file_id, var_id, i
  
  
      sizes = SHAPE(larr)
  
      ! Create the file, overwriting if it exists
      ierr = nf90_create(filename, NF90_CLOBBER, file_id)
  
      ! I don't want to bomb if there is an error, rather return to caller
      ! This is tricky to do from another sub. so I choose this instead
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
  
      ! Now I am going to do several actions, checking and printing each one
      ! I am using a loop here, to save a bit of duplication. In higher
      ! dimensions, this would really help!
  
      DO i = 1, ndims
        ierr = nf90_def_dim(file_id, dims(i), sizes(i), dim_ids(i))
        IF (ierr /= nf90_noerr) THEN
          PRINT*, TRIM(nf90_strerror(ierr))
          RETURN
        END IF
      END DO
  
      ! Define variable type, matching our array
      ierr = nf90_def_var(file_id, "grid_data", NF90_DOUBLE, dim_ids, var_id)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
  
      ! Finish defining metadata
      ierr = nf90_enddef(file_id)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
  
  
      ! Actually write the variable
      ierr = nf90_put_var(file_id, var_id, larr)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
  
      ! Close the file
      ierr = nf90_close(file_id)
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
  
    END SUBROUTINE writer_prototype
  
  
  END MODULE write_netcdf
  
  
  
  
  PROGRAM MAIN
  
    USE ISO_FORTRAN_ENV
    USE domain_tools
    USE GAUSS_SEIDEL
    USE write_netcdf
  
  
    IMPLICIT NONE
  
    INTEGER :: n_x, n_y, ierr
    REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: E_x, E_y, phi, rho
    CHARACTER(LEN=*), PARAMETER :: problem = "double"
  
  
    n_x = 100
    n_y = 100
  
  
    ALLOCATE(E_x(1:n_x, 1:n_y))
    ALLOCATE(E_y(1:n_x, 1:n_y))
  
  
    CALL get_potential(n_x, n_y, phi, problem)
    CALL get_charge_density(n_x, n_y, problem, rho)
  
    E_x = get_electric_field_x(n_x, n_y, phi)
    E_y = get_electric_field_y(n_x, n_y, phi)
  
    CALL writer_prototype(E_x, 'E_x.nc', ierr)
    CALL writer_prototype(E_y, 'E_y.nc', ierr)
    CALL writer_prototype(phi, 'phi.nc', ierr)
    CALL writer_prototype(rho, 'rho.nc', ierr)
  
  
  
  
  END PROGRAM MAIN