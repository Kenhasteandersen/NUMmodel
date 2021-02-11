#!/bin/sh
#
# For testing using: -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan

# Compile each file:
gfortran -c -fPIC  -O3 globals.f90 -o globals.o
gfortran -c -fPIC  -O3 debug.f90 -o debug.o
gfortran -c -fPIC  -O3 spectrum.f90 -o spectrum.o
gfortran -c -fPIC  -O3 generalists.f90 -o generalists.o
gfortran -c -fPIC  -O3 generalists_csp.f90 -o generalists_csp.o
gfortran -c -fPIC  -O3 copepods.f90 -o copepods.o
gfortran -c -fPIC  -O3 NUMmodel.f90 -o NUMmodel.o
gfortran -c -fPIC  -O3 NUMmodeltest.f90 -o NUMmodeltest.o
gfortran -c -fPIC  -O3 NUMmodel_wrap_colmajor.f90 -o NUMmodel_wrap_colmajor.o
# Make executable:
gfortran -fPIC  -O3  globals.o debug.o spectrum.o generalists.o generalists_csp.o copepods.o NUMmodel.o NUMmodeltest.o -o NUMmodel
# Make library
gfortran -fPIC -shared globals.o debug.o spectrum.o generalists.o generalists_csp.o copepods.o NUMmodel.o NUMmodel_wrap_colmajor.o -o NUMmodel.so

#gfortran -shared -fPIC  -o Globals.o debug.o Spectrum.o Generalists.o NUMmodel.o -o NUMmodel.so
#gfortran -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -lSystem -fPIC  -o Globals.o debug.o Spectrum.o Generalists.o NUMmodel.o -o NUMmodel
# c-wrap in two forms
#gfortran -shared -fPIC NUMmodel.o NUMmodel_wrap.f90 -o libNUMmodel.so
#gfortran -shared -fPIC -O3 NUMmodel.o NUMmodel_wrap_colmajor.f90 -o libNUMmodel_colmajor.so
