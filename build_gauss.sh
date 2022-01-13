#!/bin/bash

echo "Compiling"
#List here all the f90 files to compile,
myf90files="command_line.f90 shared_data.f90 create_axis.f90 " 
myf90files+="io.f90 velocity_verlet.f90 gauss_seidel.f90 main.f90"

### command_line.f90 not my own work, written by Chris Brady and Heather Ratcliffe, University of Warwick

#Name of compiled file
exename="laplace_solver"

# Command line variables

n_x=100
n_y=100
problem="double"


#Name of compiler
fc=gfortran

#Use nf-config to grab the compile and link flags. Backticks run command and grab output
fflags=`nf-config --fflags`
flibs=`nf-config --flibs`


#Actual compile line
$fc -std=f2008 -Wall  $fflags $myf90files $flibs -o $exename


#Report on success
if [ $? == 0 ] #if sucessful execute ising
then
  echo "Produced executable "$exename

  # runs executable using set parameters
  ./$exename n_x=$n_x n_y=$n_y problem=$problem

  echo "Python"
  python3 python_plotter.py
   
else # there is an error 
  echo "ERROR: Failed to link executable "$exename
fi

