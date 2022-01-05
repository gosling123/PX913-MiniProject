

MODULE write_netcdf

  USE ISO_FORTRAN_ENV
  USE netcdf

  IMPLICIT NONE

  CONTAINS

  ! subroutine to perform error status check
  SUBROUTINE check(status)
    INTEGER, INTENT (IN) :: status
    
    IF(status /= nf90_noerr) THEN 
      PRINT *, trim(nf90_strerror(status))
      STOP "Stopped"
    END IF
  end subroutine check  


  ! I wrote this based on the example at
  ! https://people.sc.fsu.edu/~jburkardt/f_src/netcdf/netcdf.html
  ! to get the hang of writing a simple file
  ! This page is available via the Wayback machine at
  ! https://web.archive.org/web/20190623025346/http://people.sc.fsu.edu/~jburkardt/f_src/netcdf/netcdf.html (accessedNov 2021)

  SUBROUTINE writer_prototype(filename, rho, phi, E_x, E_y, x, y, v_x, v_y, a_x, a_y, n_x, n_y, problem, ver_iter, dx, dy, d_t)

      REAL(KIND=REAL64), INTENT(IN), DIMENSION(:,:) :: rho, phi, E_x, E_y ! Electric 2D array quantaties
      REAL(KIND=REAL64), INTENT(IN), DIMENSION(:) :: x, y, v_x, v_y, a_x, a_y ! Verlet 1D array quantaties
      INTEGER, PARAMETER :: ndims = 2
      INTEGER, DIMENSION(ndims) :: size_rho, size_phi, size_E_x, size_E_y 
      INTEGER :: size_x, size_y
      INTEGER :: size_v_x, size_v_y
      INTEGER :: size_a_x, size_a_y
      INTEGER, DIMENSION(ndims) :: rho_dim_ids, phi_dim_ids, E_x_dim_ids, E_y_dim_ids
      INTEGER :: x_dim_ids, y_dim_ids
      INTEGER :: v_x_dim_ids, v_y_dim_ids
      INTEGER :: a_x_dim_ids, a_y_dim_ids
      CHARACTER(LEN=*), DIMENSION(ndims), PARAMETER :: rho_dims=(/"rho_x", "rho_y"/)
      CHARACTER(LEN=*), DIMENSION(ndims), PARAMETER :: phi_dims=(/"phi_x", "phi_y"/)
      CHARACTER(LEN=*), DIMENSION(ndims), PARAMETER :: E_x_dims=(/"E_x_x", "E_x_y"/)
      CHARACTER(LEN=*), DIMENSION(ndims), PARAMETER :: E_y_dims=(/"E_y_x", "E_y_y"/)
      CHARACTER(LEN=*), PARAMETER :: x_dims="x", y_dims="y"
      CHARACTER(LEN=*), PARAMETER :: v_x_dims="v_x", v_y_dims="v_y"
      CHARACTER(LEN=*), PARAMETER :: a_x_dims="a_x", a_y_dims="a_y"
      CHARACTER(LEN=*), INTENT(IN) :: filename
      INTEGER :: rho_var_id, phi_var_id, E_x_var_id, E_y_var_id
      INTEGER :: x_var_id, y_var_id, v_x_var_id, v_y_var_id, a_x_var_id, a_y_var_id
      INTEGER, INTENT(IN) :: n_x, n_y, ver_iter
      REAL(KIND=REAL64), INTENT(IN) :: d_t, dx, dy
      INTEGER :: k, file_id
      CHARACTER(LEN=*) :: problem
  
  
      ! define size of arrays
      size_rho = SHAPE(rho)
      size_phi = SHAPE(phi)
      size_E_x = SHAPE(E_x)
      size_E_y = SHAPE(E_y)
      size_x = SIZE(x)
      size_y = SIZE(y)
      size_v_x = SIZE(v_x)
      size_v_y = SIZE(v_y)
      size_a_x = SIZE(a_x)
      size_a_y = SIZE(a_y)
      



  
      ! Create the file, overwriting if it exists
      CALL check(nf90_create(filename, NF90_CLOBBER, file_id))
  
      ! write in the dimensions for the variables
      DO k = 1, ndims
          
        CALL check(nf90_def_dim(file_id, rho_dims(k), size_rho(k), rho_dim_ids(k)))
        CALL check(nf90_def_dim(file_id, phi_dims(k), size_phi(k), phi_dim_ids(k)))
        CALL check(nf90_def_dim(file_id, E_x_dims(k), size_E_x(k), E_x_dim_ids(k)))
        CALL check(nf90_def_dim(file_id, E_y_dims(k), size_E_y(k), E_y_dim_ids(k)))
      
      END DO


      CALL check(nf90_def_dim(file_id, x_dims, size_x, x_dim_ids))
      CALL check(nf90_def_dim(file_id, y_dims, size_y, y_dim_ids))
      CALL check(nf90_def_dim(file_id, v_x_dims, size_v_x, v_x_dim_ids))
      CALL check(nf90_def_dim(file_id, v_y_dims, size_v_y, v_y_dim_ids))
      CALL check(nf90_def_dim(file_id, a_x_dims, size_a_x, a_x_dim_ids))
      CALL check(nf90_def_dim(file_id, a_y_dims, size_a_y, a_y_dim_ids))
   
 


  
      ! Define variable type, matching the arrays
      CALL check (nf90_def_var(file_id, "rho", NF90_DOUBLE, rho_dim_ids, rho_var_id))
      CALL check (nf90_def_var(file_id, "phi", NF90_DOUBLE, phi_dim_ids, phi_var_id))
      CALL check (nf90_def_var(file_id, "E_x", NF90_DOUBLE, E_x_dim_ids, E_x_var_id))
      CALL check (nf90_def_var(file_id, "E_y", NF90_DOUBLE, E_y_dim_ids, E_y_var_id))
      CALL check (nf90_def_var(file_id, "x", NF90_DOUBLE, x_dim_ids, x_var_id))
      CALL check (nf90_def_var(file_id, "y", NF90_DOUBLE, y_dim_ids, y_var_id))
      CALL check (nf90_def_var(file_id, "v_x", NF90_DOUBLE, v_x_dim_ids, v_x_var_id))
      CALL check (nf90_def_var(file_id, "v_y", NF90_DOUBLE, v_y_dim_ids, v_y_var_id))
      CALL check (nf90_def_var(file_id, "a_x", NF90_DOUBLE, a_x_dim_ids, a_x_var_id))
      CALL check (nf90_def_var(file_id, "a_y", NF90_DOUBLE, a_y_dim_ids, a_y_var_id))
     




      CALL check(nf90_put_att(file_id, NF90_GLOBAL, "n_x", n_x))
      CALL check(nf90_put_att(file_id, NF90_GLOBAL, "n_y", n_y))
      CALL check(nf90_put_att(file_id, NF90_GLOBAL, "Prob", problem))
      CALL check(nf90_put_att(file_id, NF90_GLOBAL, "Verlet_Iteration", ver_iter))
      CALL check(nf90_put_att(file_id, NF90_GLOBAL, "Time_step", d_t))
      CALL check(nf90_put_att(file_id, NF90_GLOBAL, "Delta_x", dx))
      CALL check(nf90_put_att(file_id, NF90_GLOBAL, "Delta_y", dy))

      ! Finish defining metadata
      CALL check(nf90_enddef(file_id))
     
  
  
      ! Actually write the variables
      CALL check(nf90_put_var(file_id, rho_var_id, rho))
      CALL check(nf90_put_var(file_id, phi_var_id, phi))
      CALL check(nf90_put_var(file_id, E_x_var_id, E_x))
      CALL check(nf90_put_var(file_id, E_y_var_id, E_y))
      CALL check(nf90_put_var(file_id, x_var_id, x))
      CALL check(nf90_put_var(file_id, y_var_id, y))
      CALL check(nf90_put_var(file_id, v_x_var_id, v_x))
      CALL check(nf90_put_var(file_id, v_y_var_id, v_y))
      CALL check(nf90_put_var(file_id, a_x_var_id, a_x))
      CALL check(nf90_put_var(file_id, a_y_var_id, a_y))
      
  
      ! Close the file
      CALL check(nf90_close(file_id))

      PRINT*, "SUCESSFULLY WRITTEN NETCDF FILE"
     
  
  END SUBROUTINE writer_prototype
  
  


END MODULE write_netcdf
