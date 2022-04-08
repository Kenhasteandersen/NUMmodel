rm ..\lib\libNUMmodel_matlab.dll
make clean
rm *.mod
gfortran -c -fPIC -O2 input.f90
gfortran -c -fPIC -O2 globals.f90
gfortran -c -fPIC -O2 spectrum.f90
gfortran -c -fPIC -O2 POM.f90
gfortran -c -fPIC -O2 *.f90

gfortran -shared -fPIC -O2 *o -o libNUMmodel_matlab.dll
mv libNUMmodel_matlab.dll ../lib/.
