#
# Define lib-names depending on the system
#
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
 libMatlab = ../lib/libNUMmodel_matlab.so
 libR = ../lib/libNUMmodel_R.so
endif

ifeq ($(UNAME), Linux)
 libMatlab = ../lib/libNUMmodel_matlab.so
 libR = ../lib/libNUMmodel_R.so
endif

ifeq ($(UNAME), Windows)
 libMatlab = ../lib/libNUMmodel_matlab.dll
 libR = ../lib/libNUMmodel_R.dll
endif


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

OBJ = globals.o spectrum.o generalists.o diatoms_simple.o diatoms.o generalists_csp.o copepods.o NUMmodel.o

all: NUMmodel libMatlab libR

debugging: FLAGS=-g -check all -fpe0 -warn -traceback -debug extended

debugging:  all

%.o: %.f90
	$(FORT) -c -fPIC $(FLAGS) -o $@ $<

NUMmodel: $(OBJ)  NUMmodeltest.o
	$(FORT) $(LINK) $(FLAGS) -o $@ $^

libMatlab: $(OBJ) NUMmodel_wrap_colmajor.o
	$(FORT) -shared $(LINK) $(FLAGS) -o $(libMatlab) $^

libR: $(OBJ) NUMmodel_wrap_R.o
	$(FORT) -shared $(LINK) $(FLAGS) -o $(libR) $^

clean:
	rm -f *.o NUMmodel $(libMatlab) $(libR)
