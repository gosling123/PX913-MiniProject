MODULE VELOCITY_VERLET !Module to handle claculation of phi, E etc.

    USE ISO_FORTRAN_ENV
    USE GAUSS_SEIDEL
    use shared_data
    IMPLICIT NONE
  
    
    CONTAINS
      
    !subroutine to validate timestep, allocate arrays and set up intial condition
    SUBROUTINE setup_verlet(problem, maxsteps, stepsize)
      CHARACTER(LEN=*), INTENT(IN) :: problem !defines the set up of problem
      integer, intent(in) :: maxsteps !max number of steps to run the simulation for
      real(kind=real64) :: stepsize
    

      dt = stepsize
      max_steps = maxsteps

      allocate(r(0:maxsteps, 2),v(0:maxsteps, 2), a(0:maxsteps, 2))
      r=0.0_REAL64; v=0.0_REAL64; a=0.0_REAL64


      !Set initial conditions based on problem
      select case (problem)

        case ("null")
        r(0, :) = (/0.0_REAL64, 0.0_REAL64/)
        v(0, :) = (/0.1_REAL64, 0.1_REAL64/)
  
        case("single")  
        r(0, :) = (/0.1_REAL64, 0.0_REAL64/)
        v(0, :) = (/0.0_REAL64, 0.0_REAL64/)
  
        case("double")
        r(0, :) = (/0.0_REAL64, 0.5_REAL64/)
        v(0, :) = (/0.0_REAL64, 0.0_REAL64/)
  
        case default
            print *, "Error in setup_verlet:"
            stop "PLEASE INPUT PROBLEM AS 'null', 'single' or 'double'"

      end select
  
    END SUBROUTINE 
  
  
    !subroutine to perform the velocity verlet algorithm
    SUBROUTINE iterate_verlet(start, nsteps)
       
        INTEGER :: start, nsteps !set number of iterations
        INTEGER :: cell(2), k 

        if (.not. allocated(r)) stop "Positions not allocated for iterate_verlet, run setup_verlet first"

        if (.not. allocated(v)) stop "Velocities not allocated for iterate_verlet, run setup_verlet first"

        if (.not. allocated(a)) stop "Acceleration not allocated for iterate_verlet, run setup_verlet first"

        if (.not. allocated(E)) stop "Electical potenial not allocated for iterate_verlet, run setup_sys first"

        if (dt .eq. -1.0_real64) stop "Error in iterate_verlet, dt has not been declared"

        if (d_x .eq. -1.0_real64) stop "Error in iterate_verlet, dx has not been declared"

        if (d_y .eq. -1.0_real64) stop "Error in iterate_verlet, dy has not been declared"

        if (start < 0) stop "Error in iterate_verlet: Trying to start at negative time"

        if (start+nsteps > max_steps) stop "Error in iterate_verlet: Trying to iterate past max time" 
  
  
        
        !find indicies on the 2d grid
        cell = FLOOR((r(0, :)+1.0_REAL64)/[d_x, d_y] + 1)
        
        !initial acceleration
        a(0, :) = E(cell(1), cell(2), :)
  
    
        DO k = 1, nsteps
  
          !verlet update on posistions
          r(k, :) = r(k-1, :) + v(k-1, :)*dt + 0.50_REAL64*a(k-1, :)*dt**2.0_REAL64
  
          !find grid indicies
          cell = FLOOR((r(k, :)+1.0_REAL64)/[d_x, d_y] + 1)
          
  
          ! exit loop if out of bounds
          IF (cell(1) > n_x .or. cell(1) < 1 .or. cell(2) > n_y .or. cell(2) < 1) THEN
            print *, "particle exited at timestep", k
            tstep_end = k
            return

          end if
  
          !verlet update on accelerations
          a(k, :) = E(cell(1), cell(2), :)
  
          !verlet update on velocities
          v(k, :) = v(k-1, :) + 0.5_REAL64*dt*(a(k, :) + a(k-1, :))
  
  
        END DO
  
  
  
      END SUBROUTINE 
  
  
  
  
  
  
  
    
  END MODULE VELOCITY_VERLET
  
    
  
  
  
  