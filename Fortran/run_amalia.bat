del ..\lib\libNUMmodel_matlab.dll
del *.mod
gfortran -c -fPIC input.f90
gfortran -c -fPIC globals.f90
gfortran -c -fPIC spectrum.f90
gfortran -c -fPIC POM.f90
gfortran -c -fPIC *f90

gfortran -shared *o -o libNUMmodel_matlab.dll
move libNUMmodel_matlab.dll ../lib/.
