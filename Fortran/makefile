#
# iFortran on macos:
#
FORT=ifort
#FLAGS=-g -check all -fpe0 -warn -traceback -debug extended
FLAGS=-O2
LINK=-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
#
# gfortran on macos:
#
FORT=gfortran
#FLAGS=-g -check all -fpe0 -warn -traceback -debug extended
FLAGS=-O2
LINK=



OBJ = globals.o debug.o spectrum.o generalists.o generalists_csp.o copepods.o NUMmodel.o

%.o: %.f90
	$(FORT) -c -fPIC $(FLAGS) -o $@ $< 

NUMmodel: $(OBJ)  NUMmodeltest.o
	$(FORT) $(LINK) $(FLAGS) -o $@ $^ 

lib: $(OBJ) NUMmodel_wrap_colmajor.o
	$(FORT) -shared $(LINK) $(FLAGS) -o NUMmodel.so $^ 

clean:
	rm *.o NUMmodel NUMmodel.so


