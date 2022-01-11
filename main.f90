program main
    
    use shared_data
    use VELOCITY_VERLET
    use GAUSS_SEIDEL
    use io

    implicit NONE

    character(len = 20) :: problem = "single"


    n_x = 10
    n_y = 10

    print *, "Program Started"

    call setup_sys()
    print *, "System Set up"

    call get_charge_density(problem)
    print *, "Rho Calculated"

    call get_potential()
    print *, "Potential Calculated"

    call get_electric_field()
    print *, "Electric Field Calculated"

    call setup_verlet(problem, 1000, 0.01d0)

    print *, max_steps

    call iterate_verlet(0, 1000)

    call writer_prototype("data.ncf")

end program