FORT=ifort
GFLAGS=-g -check all -fpe0 -warn -traceback -debug extended
OFLAGS=-O3
LINK=-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
OBJ = globals.o debug.o spectrum.o generalists.o generalists_csp.o copepods.o NUMmodel.o NUMmodeltest.o

%.o: %.f90
	$(FORT) -c -fPIC -o $@ $< $(GFLAGS)

NUMmodel: $(OBJ)
	$(FORT) $(LINK) -o $@ $^ $(CFLAGS)

NUMmodel.so: $(OBJ)
	$(FORT) -shared $(LINK) -o $@ $^ $(CFLAGS)

clean:
	rm *.o NUMmodel NUMmodel.so


