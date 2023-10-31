cmake_minimum_required(VERSION 3.20)

##compiler
#set(CMAKE_C_COMPILER gcc)
set(CMAKE_Fortran_COMPILER gfortran)
#set(CMAKE_CXX_COMPILER g++)

enable_language(Fortran)
project(ref_sol Fortran)

add_library(refsol2Dhel
  ref_sol_2d_hel.f90
)

# option
target_compile_options(refsol2Dhel PRIVATE -cpp -march=native -O3 -fmax-errors=3 -ffree-line-length-none)
