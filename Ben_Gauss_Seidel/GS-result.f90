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
  
  ! allocate without ghost cells
  ALLOCATE(E_x(1:n_x, 1:n_y))
  ALLOCATE(E_y(1:n_x, 1:n_y))
  
  ! output charge density (rho) and potential (phi)
  CALL get_potential(n_x, n_y, phi, problem)
  CALL get_charge_density(n_x, n_y, problem, rho)
   
  ! use functions to output E
  E_x = get_electric_field_x(n_x, n_y, phi)
  E_y = get_electric_field_y(n_x, n_y, phi)
  
  ! Write to netcdf to plot and inspect result
  CALL writer_prototype(E_x, 'E_x.nc', ierr)
  CALL writer_prototype(E_y, 'E_y.nc', ierr)
  CALL writer_prototype(phi, 'phi.nc', ierr)
  CALL writer_prototype(rho, 'rho.nc', ierr)
  
  
  
  
END PROGRAM MAIN