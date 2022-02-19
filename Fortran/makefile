#
# Define lib-names depending on the system
#
ifeq '$(findstring ;,$(PATH))' ';'
    UNAME := Windows
else
    UNAME := $(shell uname 2>/dev/null || echo Unknown)
endif

ifeq ($(UNAME), Darwin)
 libMatlab = ../lib/libNUMmodel_matlab.dylib
 libR = ../lib/libNUMmodel_R.dylib
endif

ifeq ($(UNAME), Linux)
 libMatlab = ../lib/libNUMmodel_linux_matlab.so
 libR = ../lib/libNUMmodel_linux_R.so
endif

ifeq ($(UNAME), Windows)
 libMatlab = ../lib/libNUMmodel_matlab.dll
 libR = ../lib/libNUMmodel_R.dll
endif

#
# Compile from OSX -> windows
#
#FORT = /usr/local/Cellar/mingw-w64/9.0.0_2/toolchain-x86_64/bin/x86_64-w64-mingw32-gfortran
#libMatlab = ../lib/libNUMmodel_matlab.dll
#libR = ../lib/libNUMmodel_R.dll
#
# iFortran:
#
#FORT=ifort
#FLAGS=-O2
#LINK=-l
#
# gfortran:
#
FORT=gfortran
#FLAGS=-g3
FLAGS=-O2
LINK=
#
# nag fortran
#
#FORT = nagfor



OBJ = globals.o spectrum.o generalists.o diatoms_simple.o diatoms.o generalists_csp.o copepods.o NUMmodel.o POM.o

all: NUMmodel libMatlab libR

debugging: FLAGS=-g -check all -fpe0 -warn -traceback -debug extended

debugging:  all

%.o: %.f90
	$(FORT) -c  $(FLAGS) -o $@ $<

NUMmodel: $(OBJ)  NUMmodeltest.o
	$(FORT) $(LINK) -fPIC $(FLAGS) -o $@ $^

libMatlab: $(OBJ) NUMmodel_wrap_colmajor.o
	$(FORT) -shared $(LINK) $(FLAGS) -o $(libMatlab) $^

libR: $(OBJ) NUMmodel_wrap_R.o
	$(FORT) -shared $(LINK) $(FLAGS) -o $(libR) $^

clean:
	rm -f *.o NUMmodel $(libMatlab) $(libR)
