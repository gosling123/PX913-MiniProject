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