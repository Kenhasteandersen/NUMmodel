del ..\lib\libNUMmodel_matlab.dll
del *.mod
gfortran -c -fPIC -O2 input.f90
gfortran -c -fPIC -O2 globals.f90
gfortran -c -fPIC -O2 spectrum.f90
gfortran -c -fPIC -O2 POM.f90
gfortran -c -fPIC -O2 *f90

gfortran -shared *o -o libNUMmodel_matlab.dll
move libNUMmodel_matlab.dll ../lib/.
