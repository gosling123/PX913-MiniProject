module io

    USE ISO_FORTRAN_ENV
    USE netcdf
    use shared_data
    implicit NONE

    CONTAINS
  
    !************************************************************************
    !> setup_sys
    !!
    !! subroutine to perform error status check
    !!
    !!@param status : status code for error
    !************************************************************************
    !
    SUBROUTINE check(status)
      INTEGER, INTENT (IN) :: status
      
      IF(status /= nf90_noerr) THEN 
        PRINT *, trim(nf90_strerror(status))
        STOP ""
      END IF
    end subroutine check  
    
    !************************************************************************
    !> writer prototype
    !! subroutine to perform error status check
    !!
    !! I wrote this based on the example at
    !! https://people.sc.fsu.edu/~jburkardt/f_src/netcdf/netcdf.html
    !! to get the hang of writing a simple file
    !! This page is available via the Wayback machine at
    !! https://web.archive.org/web/20190623025346/http://people.sc.fsu.edu/~jburkardt/f_src/netcdf/netcdf.html (accessedNov 2021)
    !!
    !!@param filename : filename to store data in
    !!@param problem : initial problem setup
    !************************************************************************
    !
    SUBROUTINE writer_prototype(filename, problem)
  
        INTEGER, DIMENSION(2) :: size_rho, size_phi, size_r, size_v, size_a
        INTEGER, DIMENSION(2) :: rho_dim_ids, phi_dim_ids, r_dim_ids, v_dim_ids, a_dim_ids
        INTEGER, DIMENSION(3) :: size_E
        INTEGER, DIMENSION(3) :: E_dim_ids
        CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: rho_dims=(/"rho_x", "rho_y"/)
        CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: phi_dims=(/"phi_x", "phi_y"/)
        CHARACTER(LEN=*), DIMENSION(3), PARAMETER :: E_dims=(/"E_x", "E_y", "d_E"/)
        CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: r_dims=(/"r_x", "r_y"/)
        CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: v_dims=(/"v_x", "v_y"/)
        CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: a_dims=(/"a_x", "a_y"/)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER :: rho_var_id, phi_var_id, E_var_id
        INTEGER :: r_var_id, v_var_id,  a_var_id
        INTEGER :: k, file_id
        CHARACTER(LEN=*) :: problem 
    
    
        ! define size of arrays
        size_rho = SHAPE(rho(1:n_x, 1:n_y))
        size_phi = SHAPE(phi(1:n_x, 1:n_y))
        size_E = SHAPE(E)
        size_r = SHAPE(r)
        size_v = SHAPE(v)
        size_a = SHAPE(a)
   
        print *, size_r
        ! Create the file, overwriting if it exists
        CALL check(nf90_create(filename, NF90_CLOBBER, file_id))
    
        ! write in the dimensions for the variables
        DO k = 1, 2

          CALL check(nf90_def_dim(file_id, rho_dims(k), size_rho(k), rho_dim_ids(k)))
          CALL check(nf90_def_dim(file_id, phi_dims(k), size_phi(k), phi_dim_ids(k)))
          CALL check(nf90_def_dim(file_id, E_dims(k), size_E(k), E_dim_ids(k)))
          CALL check(nf90_def_dim(file_id, r_dims(k), size_r(k), r_dim_ids(k)))
          CALL check(nf90_def_dim(file_id, v_dims(k), size_v(k), v_dim_ids(k)))
          CALL check(nf90_def_dim(file_id, a_dims(k), size_a(k), a_dim_ids(k)))
        
        END DO
  
  
        CALL check(nf90_def_dim(file_id, E_dims(3), size_E(3), E_dim_ids(3)))

        print *, "Dimensions Written"

    
        ! Define variable type, matching the arrays
        CALL check (nf90_def_var(file_id, "rho", NF90_DOUBLE, rho_dim_ids, rho_var_id))
        CALL check (nf90_def_var(file_id, "phi", NF90_DOUBLE, phi_dim_ids, phi_var_id))
        CALL check (nf90_def_var(file_id, "E", NF90_DOUBLE, E_dim_ids, E_var_id))
        CALL check (nf90_def_var(file_id, "r", NF90_DOUBLE, r_dim_ids, r_var_id))
        CALL check (nf90_def_var(file_id, "v", NF90_DOUBLE, v_dim_ids, v_var_id))
        CALL check (nf90_def_var(file_id, "a", NF90_DOUBLE, a_dim_ids, a_var_id))
       
        print *, "Variable types defined"
  
  
        !Setting up meta variables of simulation
        CALL check(nf90_put_att(file_id, NF90_GLOBAL, "Prob", problem))
        CALL check(nf90_put_att(file_id, NF90_GLOBAL, "Verlet_Iteration", viter))
        CALL check(nf90_put_att(file_id, NF90_GLOBAL, "Delta_x", d_x))
        CALL check(nf90_put_att(file_id, NF90_GLOBAL, "Delta_y", d_y))

        ! Finish defining metadata
        CALL check(nf90_enddef(file_id))
       
    
    
        ! Actually write the variables
        CALL check(nf90_put_var(file_id, rho_var_id, rho(1:n_x, 1:n_y)))
        CALL check(nf90_put_var(file_id, phi_var_id, phi(1:n_x, 1:n_y)))
        CALL check(nf90_put_var(file_id, E_var_id, E))
        CALL check(nf90_put_var(file_id, r_var_id, r))
        CALL check(nf90_put_var(file_id, v_var_id, v))
        CALL check(nf90_put_var(file_id, a_var_id, a))
        
    
        ! Close the file
        CALL check(nf90_close(file_id))
  
        PRINT*, "SUCESSFULLY WRITTEN NETCDF FILE"
       
    
    END SUBROUTINE writer_prototype
    
end module  
