MODULE VELOCITY_VERLET !Module to handle claculation of phi, E etc.

    USE ISO_FORTRAN_ENV
    USE GAUSS_SEIDEL
    use shared_data
    IMPLICIT NONE
  
    
    CONTAINS
    
    !************************************************************************
    !> setup_verlet
    !!
    !!subroutine to validate timestep, allocate arrays and set up intial condition
    !!
    !! @param problem : initial problem setup
    !! @param maxsteps : maximum number of steps for simulation
    !! @param stepsize : time interval between two iterations (dt)
    !************************************************************************
    !    
    SUBROUTINE setup_verlet(problem, maxsteps, stepsize)
      CHARACTER(LEN=*), INTENT(IN) :: problem 
      integer, intent(in) :: maxsteps
      real(kind=real64) :: stepsize
    

      d_t = stepsize
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
  
  
    !************************************************************************
    !> iterative_verlet
    !!
    !!subroutine to perform the velocity verlet algorithm
    !!
    !! @param start : start index for iterating verlet
    !! @param nsteps : number of steps for simulation after start
    !************************************************************************
    !
    SUBROUTINE iterate_verlet(start, nsteps)
       
        INTEGER :: start, nsteps 
        INTEGER :: cell(2), k 

        if (.not. allocated(r)) stop "Positions not allocated for iterate_verlet, run setup_verlet first"

        if (.not. allocated(v)) stop "Velocities not allocated for iterate_verlet, run setup_verlet first"

        if (.not. allocated(a)) stop "Acceleration not allocated for iterate_verlet, run setup_verlet first"

        if (.not. allocated(E)) stop "Electical potenial not allocated for iterate_verlet, run setup_sys first"

        if (d_t .eq. -1.0_real64) stop "Error in iterate_verlet, d_t has not been declared"

        if (start < 0) stop "Error in iterate_verlet: Trying to start at negative time"

        if (start+nsteps > max_steps) stop "Error in iterate_verlet: Trying to iterate past max time" 
  
  
        
        !find indicies on the 2d grid
        cell = FLOOR(  (r(0, :)+1.0_REAL64)  /[d_x, d_y] + 1)
      
        !initial acceleration
        a(0, :) = q/m *E(cell(1), cell(2), :)
  
        DO k = 1, nsteps
  
          !verlet update on positions
          r(k, :) = r(k-1, :) + v(k-1, :)*d_t + 0.50_REAL64*a(k-1, :)*d_t**2.0_REAL64
          
          !find grid indicies
          cell = FLOOR((r(k, :)+1.0_REAL64)/[d_x, d_y] + 1)
          
          ! exit loop if out of bounds
          IF (cell(1) > n_x .or. cell(1) < 1 .or. cell(2) > n_y .or. cell(2) < 1) THEN
            print *, "particle exited at timestep", k
            viter = k
            return

          end if
  
          !verlet update on accelerations
          a(k, :) = q/m*E(cell(1), cell(2), :)
    
          !verlet update on velocities
          v(k, :) = v(k-1, :) + 0.5_REAL64*d_t*(a(k, :) + a(k-1, :))

        END DO
  
        viter = k
  
      END SUBROUTINE
    
  END MODULE VELOCITY_VERLET
  
    
  
  
  
  