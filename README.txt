READ ME


This is a short text describing the functions of the different source files


- main.f90: Main program, calls on functions from all other .f90 files

- shared_data.f90: Contains the shared_data module for all other .f90 modules

- gauss_seidel.f90: Contains subroutines to set up system and preform gauss seidel

- velocity_verlet.f90: Contains subroutines to set up and execute velocity verlet integrator

- io.f90: Writes a NetCDF file wth trajectory and system data

- create_axis.f90: File given by Christopher Brady and Heather Ratcliffe for setting up grid arrays

- python_plotter.py: Reads NetCDF file and outputs visualization 

- build_gauss.sh : Compiles and Runs main.f90 with its used modules. Then it runs the python_plotter.py
