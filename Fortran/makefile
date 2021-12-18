#INC = /Library/Frameworks/R.framework/Resources # For C-wrapper for R

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
# Link needed on Macos:
#LINK=-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
#LINK=
#
# lfortran
#

OBJ = globals.o debug.o spectrum.o generalists.o generalists_csp.o diatoms.o diatoms_simple.o copepods.o NUMmodel.o

all: NUMmodel libMatlab libR

#debugging: FLAGS = -g3 -fbounds-check
debugging: FLAGS=-g -check all -fpe0 -warn -traceback -debug extended

debugging: clean all

%.o: %.f90
	$(FORT) -c -fPIC $(FLAGS) -o $@ $<

NUMmodel: $(OBJ)  NUMmodeltest.o
	$(FORT) $(LINK) $(FLAGS) -o $@ $^

libMatlab: $(OBJ) NUMmodel_wrap_colmajor.o
	$(FORT) -shared $(LINK) $(FLAGS) -o NUMmodel_matlab.so $^

libR: $(OBJ) NUMmodel_wrap_R.o
	$(FORT) -shared $(LINK) $(FLAGS) -o NUMmodel_R.so $^

clean:
	rm -f *.o NUMmodel NUMmodel_matlab.so NUMmodel_R.so
