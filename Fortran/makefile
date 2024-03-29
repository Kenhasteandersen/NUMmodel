#
# Try to determine the operating system:
#
ifeq '$(findstring ;,$(PATH))' ';'
    UNAME := Windows
else
    UNAME := $(shell uname 2>/dev/null || echo Unknown)
    ifeq ($(UNAME), Unknown) 
        UNAME := Windows # Assuming windows if we cannot find a uname
        $(warning    Cannot determine system; assuming Windows)
    endif
endif

$(info    System is $(UNAME))

# Check for other windows name
ifeq ($(findstring NT-,$(UNAME)),NT-)
    # Found
    UNAME=Windows
endif
#
# Define lib-names depending on the system:
#
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
# iFortran:
#
#FORT=ifort
#FLAGS=-fPIC#-O2
#LINK= #-l
#
# gfortran:
#
FORT=gfortran
#FLAGS=-g3
FLAGS= -fPIC -O2
#LINK=
#
# nag fortran
#
#FORT = nagfor



OBJ = input.o globals.o spectrum.o generalists.o diatoms_simple.o diatoms.o generalists_csp.o copepods.o POM.o NUMmodel.o 

all: NUMmodel libMatlab libR

debugging: FLAGS=-fPIC -g3 -fbacktrace#-debug extended # -check all -fpe0 -warn -traceback 

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
