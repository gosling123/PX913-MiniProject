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

  SUBROUTINE check_error(ierr)
    INTEGER,INTENT(IN) :: ierr

    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      STOP 'Stopped'
    END IF

  END SUBROUTINE check_error
  ! The below subroutine takes a 2-D array larr a 1-D array time_arr
  ! a rundata(for global metadata) type and writes it into a Netcdf file 
  ! @param larr 2-D array of type integer
  ! @param time_arr 1-D array of type real
  ! @param run_data rundata type structure that contains information about the code
  ! @param filename name of the file to create and store data in 
  ! @param ierr to store return error code
  SUBROUTINE writer_prototype(filename,rho,phi,E_x,E_y,x,y,v_x,v_y,a_x,a_y,n_x,n_y,problem,ver_iter,dx,dy,dt)

    REAL(REAL64),INTENT(IN),DIMENSION(:,:)  :: rho,phi,E_x,E_y
    REAL(REAL64),INTENT(IN),DIMENSION(:)    :: x,y,v_x,v_y,a_x,a_y
    INTEGER, PARAMETER                      :: ndims = 2
    INTEGER, DIMENSION(ndims)               :: rho_size,phi_size,E_x_size,E_y_size
    INTEGER, DIMENSION(ndims)               :: rho_dim_ids,phi_dim_ids,E_x_dim_ids,E_y_dim_ids
    INTEGER                                 :: x_size,y_size,v_x_size,v_y_size,a_x_size,a_y_size 
    CHARACTER(LEN=*), DIMENSION(ndims),PARAMETER      :: rho_dims=(/"rho_x", "rho_y"/),phi_dims=(/"phi_x","phi_y"/)
    CHARACTER(LEN=*), DIMENSION(ndims),PARAMETER      :: E_x_dims=(/"E_x_x", "E_x_y"/),E_y_dims=(/"E_y_x","E_y_y"/)
    CHARACTER(LEN=*), PARAMETER             :: x_dim = "x", y_dim="y", v_x_dim="v_x", v_y_dim="v_y", a_x_dim="a_x", a_y_dim="a_y"
    INTEGER                                 :: x_dim_id,y_dim_id,v_x_dim_id,v_y_dim_id,a_x_dim_id,a_y_dim_id
    CHARACTER(LEN=*), INTENT(IN)            :: filename
    INTEGER                                 :: file_id
    INTEGER                                 :: rho_id,phi_id,E_x_id,E_y_id
    INTEGER                                 :: x_id,y_id,v_x_id,v_y_id,a_x_id,a_y_id
    INTEGER                                 :: ierr, i
    INTEGER,INTENT(IN)                      :: n_x,n_y,ver_iter
    REAL(REAL64),INTENT(IN)                 :: dx,dy,dt
    CHARACTER(LEN=*),INTENT(IN)             :: problem

    rho_size = SHAPE(rho)
    phi_size = SHAPE(phi)
    E_x_size = SHAPE(E_x)
    E_y_size = SHAPE(E_y)
    
    x_size = SIZE(x)
    y_size = SIZE(y)
    v_x_size = SIZE(v_x)
    v_y_size = SIZE(v_y)
    a_x_size = SIZE(a_x)
    a_y_size = SIZE(a_y)

    ! Create the file, overwriting if it exists
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    CALL check_error(ierr)

    !Defining global attributes
    ierr = nf90_put_att(file_id, NF90_GLOBAL, "n_x", n_x)
    CALL check_error(ierr)

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "n_y", n_y)
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

      ierr = nf90_def_dim(file_id, E_x_dims(i), E_x_size(i), E_x_dim_ids(i))
      CALL check_error(ierr)

      ierr = nf90_def_dim(file_id, E_y_dims(i), E_y_size(i), E_y_dim_ids(i))
      CALL check_error(ierr)

    END DO

    ! Defining dimension for 1D arrays
    ierr = nf90_def_dim(file_id, x_dim, x_size, x_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, y_dim, y_size, y_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, v_x_dim, v_x_size, v_x_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, v_y_dim, v_y_size, v_y_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, a_x_dim, a_x_size, a_x_dim_id)
    CALL check_error(ierr)

    ierr = nf90_def_dim(file_id, a_y_dim, a_y_size, a_y_dim_id)
    CALL check_error(ierr)

    ! Define variable type for data
    ierr = nf90_def_var(file_id, "rho", NF90_DOUBLE, rho_dim_ids, rho_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "phi", NF90_DOUBLE, phi_dim_ids, phi_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "E_x", NF90_DOUBLE, E_x_dim_ids, E_x_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "E_y", NF90_DOUBLE, E_y_dim_ids, E_y_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "x", NF90_DOUBLE, x_dim_id, x_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "y", NF90_DOUBLE, y_dim_id, y_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "v_x", NF90_DOUBLE, v_x_dim_id, v_x_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "v_y", NF90_DOUBLE, v_y_dim_id, v_y_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "a_x", NF90_DOUBLE, a_x_dim_id, a_x_id)
    CALL check_error(ierr)

    ierr = nf90_def_var(file_id, "a_y", NF90_DOUBLE, a_y_dim_id, a_y_id)
    CALL check_error(ierr)

    ! Finish defining metadata
    ierr = nf90_enddef(file_id)
    CALL check_error(ierr)

    ! Write the variables
    ierr = nf90_put_var(file_id, rho_id, rho)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, phi_id, phi)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, E_x_id, E_x)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, E_y_id, E_y)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, x_id, x)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, y_id, y)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, v_x_id, v_x)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, v_y_id, v_y)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, a_x_id, a_x)
    CALL check_error(ierr)

    ierr = nf90_put_var(file_id, a_y_id, a_y)
    CALL check_error(ierr)

    ! Close the file
    ierr = nf90_close(file_id)
    CALL check_error(ierr)

  END SUBROUTINE writer_prototype


END MODULE write_netcdf
