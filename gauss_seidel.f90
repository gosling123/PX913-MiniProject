
! Module to setup a lattice model of the system.
!
! Contains subroutines for calculating the scalar potential,
!charge density and electci field lines

MODULE GAUSS_SEIDEL 

  USE ISO_FORTRAN_ENV
  USE domain_tools
  use shared_data
  IMPLICIT NONE

  
  CONTAINS

  !Subroutine to allocate memory  and to genearte grid structure 
  ! using create axis supplied by Chris Brady and Heather 
  ! Ratcliffe (University of Warwick)

  ! System Properties Variables are all set to 0 here

  SUBROUTINE setup_sys()

    !Creates grid_X-axis
    CALL create_axis(grid_X, nels=n_x, axis_range=[-1.0_real64, 1.0_real64], nghosts=1, delta=d_x)
    print *, "X axis generated"
    

    !Creates grid_Y-axis
    CALL create_axis(grid_Y, nels=n_y, axis_range=[-1.0_real64, 1.0_real64], nghosts=1, delta=d_y)
    print *, "Y axis generated"

    
    !Allocating and setting to 0
    allocate(rho(0:n_x+1, 0:n_y+1), phi(0:n_x+1, 0:n_y+1))
    allocate(E(n_x, n_y, 2))

    rho = 0; phi = 0; E = 0;

  END SUBROUTINE setup_sys


  !Subroutine to handle the generation of selected charge densities
  SUBROUTINE get_charge_density(problem)

    INTEGER :: i, j
    CHARACTER(LEN=*), INTENT(IN) :: problem !defines the set up of rho
    REAL(KIND=REAL64) :: exponent_1, exponent_2
      

    if (.not. allocated(rho)) then
      stop "Error in get_charge_density: Charge density not allocated yet, please call setup_sys first"
    end if


    select case (problem)

      ! null -- charge density is equal to zero in all space
      case ("null") 
      
      !Rho already initialized to be 0, exit subroutine
      return


      !single -- charge density has spatial variation defined by a single Gaussian function
      case ("single") 
                
      DO j=1, n_y
        DO i=1, n_x

          exponent_1 = ((grid_X(i)/0.1)**2 + (grid_Y(j)/0.1)**2)
          rho(i,j) = EXP(-1.0_REAL64*exponent_1)
        
        END DO 
      END DO 



    !double -- charge density has spatial variation defined by a double Gaussian function  
      case ("double")

      DO j=1, n_y
        DO i=1, n_x

          exponent_1 = -1.0_REAL64*(((grid_X(i)+0.25)/0.1)**2 + ((grid_Y(j)+0.25)/0.1)**2) 
          exponent_2 = -1.0_REAL64*(((grid_X(i)-0.75)/0.2)**2 + ((grid_Y(j)-0.75)/0.2)**2)
          rho(i,j) = EXP(exponent_1) + EXP(exponent_2)
        
        END DO 
      END DO

      ! condition to enforce one of the three setups
      case default
      PRINT*, "PLEASE INPUT PROBLEM AS 'null', 'single' or 'double'"
      STOP "Stopped"

      END SELECT

  END SUBROUTINE get_charge_density


  !subroutine to calculate electric potential from Poission's Equation using Gauss-Seidel iteration
  SUBROUTINE get_potential()

    INTEGER :: i, j, k, N
    REAL(KIND=REAL64) :: x_part, y_part, denominator
    REAL(KIND=REAL64) :: e_tot, d_rms, x_FD, y_FD, x_error, y_error, d

    REAL(KIND=REAL64), PARAMETER :: tol = 0.00001_REAL64
    INTEGER, PARAMETER :: iterations = 1000000


    if (.not. allocated(phi) .or. (.not. allocated(rho))) then
      stop "Error in get_potential: phi or rho is not allocated yet, please call setup_sys first"
    end if

    !Denominator of Gauss-Seidel iteration
    denominator = 2/(d_x)**2 + 2/(d_y)**2  
     
    !e_tot - error defined as the sum of the differnces between the finite difference form of the laplacian of phi
    ! and the value of rho at each point on the grid.
    e_tot = 1.0_REAL64

    !d_rms - root mean square of the calculated error in phi
    d_rms = 1.0_REAL64

    N = n_x*n_y 


    IF (all(rho .eq. 0)) THEN
      return

    ELSE
      !Start of Gauss-Seidel process
      DO k=1, iterations

        IF (ABS(e_tot/d_rms) <= tol) then
          print *, "Gauss Seidel Converged Successfully after", k, "iterations" 
          RETURN !condition for breaking loop and outputting a converged value of phi

        end if 
        e_tot = 0.0_REAL64
        d_rms = 0.0_REAL64
  
        DO j=1, n_y
          DO i=1, n_x
            
            x_part = (phi(i+1,j) + phi(i-1,j))/(d_x)**2
            y_part = (phi(i,j+1) + phi(i,j-1))/(d_y)**2

            phi(i,j) = -1.0_REAL64*(rho(i,j)-x_part-y_part)/denominator !Gauss-Seidel Iteartion of phi

          END DO 
        END DO


        DO j=1, n_y
          DO i=1, n_x

            x_FD = (phi(i+1,j) - 2*phi(i,j) + phi(i-1,j))/(d_x)**2 !second order x deriavtive             
            y_FD = (phi(i,j+1) - 2*phi(i,j) + phi(i,j-1))/(d_y)**2 !second order y deriavtive

            x_error =  x_FD
            y_error =  y_FD

            e_tot = e_tot + ABS(x_error + y_error - rho(i,j))  ! summing to find total error

            d = SQRT(1.0_REAL64/REAL(N, REAL64))*SQRT(x_error**2 + y_error**2)

            d_rms = d_rms + d ! summing to find total root mean square

          END DO 
        END DO

        PRINT*,'Iteration:', k
        print*, 'e_tot', e_tot,'d_rms', d_rms, 'Condition', e_tot/d_rms

      END DO 

      ! Warning if phi has failed to converge
      PRINT *, "WARNING, GAUSS SEIDEL NOT CONVERGED"

    END IF

  END SUBROUTINE get_potential

  !subroutine to calculate the electric field from phi t
  subroutine get_electric_field()
    INTEGER :: i, j

    if (.not. allocated(phi) .or. (.not. allocated(rho))) then
      stop "Error in get_potential: phi or rho is not allocated yet, please call setup_sys first"
    end if

    DO j=1, n_y
      DO i=1, n_x

        E(i,j, 0) = (phi(i+1,j) - phi(i-1,j))/(2*d_x) !Electric field in x
        E(i,j, 1) = (phi(i,j+1) - phi(i,j-1))/(2*d_y) !Electric field in y
      
      END DO 
    END DO

  END subroutine

  
END MODULE GAUSS_SEIDEL

  



