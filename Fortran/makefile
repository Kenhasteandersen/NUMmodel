FORT=ifort
FLAGS=-g -check all -fpe0 -warn -traceback -debug extended
#FLAGS=-O3
LINK=-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
OBJ = globals.o debug.o spectrum.o generalists.o generalists_csp.o copepods.o NUMmodel.o

%.o: %.f90
	$(FORT) -c -fPIC -o $@ $< $(FLAGS)

NUMmodel: $(OBJ)  NUMmodeltest.o
	$(FORT) $(LINK) -o $@ $^ $(FLAGS)

lib: $(OBJ) NUMmodel_wrap_colmajor.o
	$(FORT) -shared $(LINK) -o NUMmodel.so $^ $(FLAGS)

clean:
	rm *.o NUMmodel NUMmodel.so


