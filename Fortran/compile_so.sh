#!/bin/sh
# Compiles various shared libraries.

gfortran -c -fPIC NUMmodel.f90 -o NUMmodel.o
# c-wrap in two forms
#gfortran -shared -fPIC NUMmodel.o NUMmodel_wrap.f90 -o libNUMmodel.so
gfortran -shared -fPIC NUMmodel.o NUMmodel_wrap_colmajor.f90 -o libNUMmodel_colmajor.so
# c-wrap with openMP
#gfortran -fopenmp -c -fPIC NUMmodel.f90 -o NUMmodel_openMP.o
#gfortran -fopenmp -shared -fPIC NUMmodel_openMP.o NUMmodel_wrap.f90 -o libNUMmodel_openMP.so
#gfortran -fopenmp -shared -fPIC NUMmodel_openMP.o NUMmodel_wrap_colmajor.f90 -o libNUMmodel_colmajor_openMP.so

# compile the c-runner using above shared lib
#gfortran -L. -lNUMmodel run_NUMmodel_from_c.c -o cout_so
#gfortran -L. -lNUMmodel_openMP run_NUMmodel_from_c.c -o cout_openmp_so

# compile the Fortran runner using a shared lib
#gfortran -shared -fPIC NUMmodel.f90 -o libNUMmodel_fortran.so
#gfortran -L. -lNUMmodel_fortran run_NUMmodel_from_fortran.f90 -o fout_so

# compile fortran directly to a shared lib
#gfortran -shared -fPIC NUMmodel.f90 -o libNUMmodel_f90.so
#gfortran          -shared -fPIC NUMmodel_wrap_f90.f90 NUMmodel.o        -o libNUMmodel_f90.so
#gfortran -fopenmp -shared -fPIC NUMmodel_wrap_f90.f90 NUMmodel_openMP.o -o libNUMmodel_f90_openMP.so

