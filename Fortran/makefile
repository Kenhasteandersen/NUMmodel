#INC = /Library/Frameworks/R.framework/Resources # For C-wrapper for R

#
# iFortran on macos:
#
FORT=ifort
FLAGS=-g -check all -fpe0 -warn -traceback -debug extended
#FLAGS=-O2
#LINK=
#
# gfortran on macos:
#
#FORT=gfortran
#FLAGS=-g3
#FLAGS=-O2
#LINK=-L/Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk/usr/lib
#LINK=



OBJ = globals.o debug.o spectrum.o generalists.o generalists_csp.o copepods.o NUMmodel.o

all: NUMmodel libMatlab libR

%.o: %.f90
	$(FORT) -c -fPIC $(FLAGS) -o $@ $<

NUMmodel: $(OBJ)  NUMmodeltest.o
	$(FORT) $(LINK) $(FLAGS) -o $@ $^

libMatlab: $(OBJ) NUMmodel_wrap_colmajor.o
	$(FORT) -shared $(LINK) $(FLAGS) -o NUMmodel_matlab.so $^

libR: $(OBJ) NUMmodel_wrap_R.o
	$(FORT) -shared $(LINK) $(FLAGS) -o NUMmodel_R.so $^

clean:
	rm *.o NUMmodel NUMmodel_matlab.so NUMmodel_R.so
