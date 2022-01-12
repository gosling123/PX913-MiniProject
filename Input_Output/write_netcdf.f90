! The below code builds upon the write_netcdf_array.f90
! provided in Worskshop7 and the basic structure for opening,
! writing and closing remains the same.
! This code was itself based on the example at
! https://people.sc.fsu.edu/~jburkardt/f_src/netcdf/netcdf.html

MODULE write_netcdf

  USE ISO_FORTRAN_ENV
  USE netcdf

  IMPLICIT NONE

  CONTAINS

  ! The below helper subroutine checks for error for any
  ! of the NetCDF calls we place and stops the program 
  ! in case an error does exist
  !@param ierr error code for NetCDF subroutine
  SUBROUTINE check_error(ierr)
    INTEGER,INTENT(IN) :: ierr

    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      STOP 'Stopped'
    END IF

  END SUBROUTINE check_error

  ! The below subroutine takes a 2-D array larr a 1-D array time_arr
  ! @param filename name of the file to create and store data in
  ! @param rho 2D array containing data for charge density
  ! @param phi 2D array containing data for scalar potential
  ! @param Ex 2D array containing data for x component for electric field
  ! @param Ey 2D array containing data for y component for electric field
  ! @param x 1D array containing x-axis position of particle at each timestep
  ! @param y 1D array containing y-axis position of particle at each timestep
  ! @param vx 1D array containing x-component of velocity of particle at each timestep
  ! @param vy 1D array containing y-component of velocity of particle at each timestep
  ! @param ax 1D array containing x-component of acceleration of particle at each timestep
  ! @param ay 1D array containing y-component of acceleration of particle at each timestep
  ! @param nx number of grid elements in the x-axis(global)
  ! @param ny number of grid elements in the y-axis(global)
  ! @param problem the user-defined initial conditions(global)
  ! @param ver_iter Number of iterations done for velocity verlet(global)
  ! @param dx size of grid discretization in x-direction(global)
  ! @param dy size of grid discretization in y-direction(global)
  ! @param dt Timestep used for velocity verlet(global)
  SUBROUTINE writer_prototype(filename,rho,phi,Ex,Ey,x,y,vx,vy,ax,ay,nx,ny,problem,ver_iter,dx,dy,dt)

    REAL(REAL64),INTENT(IN),DIMENSION(:,:)        :: rho,phi,Ex,Ey
    REAL(REAL64),INTENT(IN),DIMENSION(:)          :: x,y,vx,vy,ax,ay
    INTEGER, PARAMETER                            :: ndims = 2
    INTEGER, DIMENSION(ndims)                     :: rho_size,phi_size,Ex_size,Ey_size
    INTEGER, DIMENSION(ndims)                     :: rho_dim_ids,phi_dim_ids,Ex_dim_ids,Ey_dim_ids
    INTEGER                                       :: x_size,y_size,vx_size,vy_size,ax_size,ay_size 
    CHARACTER(LEN=*), DIMENSION(ndims),PARAMETER  :: rho_dims=(/"rho_x", "rho_y"/),phi_dims=(/"phi_x","phi_y"/)
    CHARACTER(LEN=*), DIMENSION(ndims),PARAMETER  :: Ex_dims=(/"Ex_x", "Ex_y"/),Ey_dims=(/"Ey_x","Ey_y"/)
    CHARACTER(LEN=*), PARAMETER                   :: x_dim = "x", y_dim="y", vx_dim="vx", vy_dim="vy", ax_dim="ax", ay_dim="ay"
    INTEGER                                       :: x_dim_id,y_dim_id,vx_dim_id,vy_dim_id,ax_dim_id,ay_dim_id
    CHARACTER(LEN=*), INTENT(IN)                  :: filename
    INTEGER                                       :: file_id
    INTEGER                                       :: rho_id,phi_id,Ex_id,Ey_id
    INTEGER                                       :: x_id,y_id,vx_id,vy_id,ax_id,ay_id
    INTEGER                                       :: ierr, i
    INTEGER,INTENT(IN)                            :: nx,ny,ver_iter
    REAL(REAL64),INTENT(IN)                       :: dx,dy,dt
    CHARACTER(LEN=*),INTENT(IN)                   :: problem

    ! Getting shapes and sizes for the variables to be written
    rho_size = SHAPE(rho)
    phi_size = SHAPE(phi)
    Ex_size = SHAPE(Ex)
    Ey_size = SHAPE(Ey)
    
    x_size = SIZE(x)
    y_size = SIZE(y)
    vx_size = SIZE(vx)
    vy_size = SIZE(vy)
    ax_size = SIZE(ax)
    ay_size = SIZE(ay)

    ! Create the file, overwriting if it exists
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    CALL check_error(ierr)

    !Defining global attributes
    ierr = nf90_put_att(file_id, NF90_GLOBAL, "nx", nx)
    CALL check_error(ierr)

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "ny", ny)
    CALL check_error(ierr)

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "Problem", problem)
    CALL check_error(ierr)

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "Verlet Integration", ver_iter)
    CALL check_error(ierr)

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "Time Step", dt)
    CALL check_error(ierr)

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "dx", dx)
    CALL check_error(ierr)
    
    ierr = nf90_put_att(file_id, NF90_GLOBAL, "dy", dy)
    CALL check_error(ierr)

    ! Defining dimensions for 2d arrays
    DO i = 1, ndims

      ierr = nf90_def_dim(file_id, rho_dims(i), rho_size(i), rho_dim_ids(i))
      CALL check_error(ierr)

      ierr = nf90_def_dim(file_id, phi_dims(i), phi_size(i), phi_dim_ids(i))
      CALL check_error(ierr)

      ierr = nf90_def_dim(file_id, Ex_dims(i), Ex_size(i), Ex_dim_ids(i))
      CALL check_error(ierr)

      ierr = nf90_def_dim(file_id, Ey_dims(i), Ey_size(i), Ey_dim_ids(i))
      CALL check_error(ierr)

    END DO

    ! Defining dimension for 1D arrays
    ierr = nf90_def_dim(file_id, x_dim, x_size, x_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, y_dim, y_size, y_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, vx_dim, vx_size, vx_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, vy_dim, vy_size, vy_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, ax_dim, ax_size, ax_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, ay_dim, ay_size, ay_dim_id)
    CALL check_error(ierr)

    ! Define variable type for data
    ierr = nf90_def_var(file_id, "rho", NF90_DOUBLE, rho_dim_ids, rho_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "phi", NF90_DOUBLE, phi_dim_ids, phi_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "Ex", NF90_DOUBLE, Ex_dim_ids, Ex_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "Ey", NF90_DOUBLE, Ey_dim_ids, Ey_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "x", NF90_DOUBLE, x_dim_id, x_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "y", NF90_DOUBLE, y_dim_id, y_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "vx", NF90_DOUBLE, vx_dim_id, vx_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "vy", NF90_DOUBLE, vy_dim_id, vy_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "ax", NF90_DOUBLE, ax_dim_id, ax_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "ay", NF90_DOUBLE, ay_dim_id, ay_id)
    CALL check_error(ierr)

    ! Finish defining metadata
    ierr = nf90_enddef(file_id)
    CALL check_error(ierr)

    ! Write the variables
    ierr = nf90_put_var(file_id, rho_id, rho)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, phi_id, phi)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, Ex_id, Ex)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, Ey_id, Ey)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, x_id, x)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, y_id, y)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, vx_id, vx)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, vy_id, vy)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, ax_id, ax)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, ay_id, ay)
    CALL check_error(ierr)

    ! Close the file
    ierr = nf90_close(file_id)
    CALL check_error(ierr)

  END SUBROUTINE writer_prototype


END MODULE write_netcdf
