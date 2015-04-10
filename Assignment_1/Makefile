FC = gfortran

FCFLAGS = -g -fbounds-check
FCFLAGS = -O2

PROGRAMS = MyExe

#all:	$(PROGRAMS)
#makedepf90
#creates a dependency to main, that main is dependent on parameters
MyExe: parameters.o main.o Ur.o f.o potential.o force.o random.o euler.o Pu.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)	

#To make an executable, one first have to make a code file .o which has
#the same dependencies as the executable test
main.o: main.f95 parameters.mod Ur.f95 f.f95 potential.f95 force.f95 random.mod euler.f95 Pu.f95
	$(FC) $(FCFLAGS) -c $<

# The .mod is compiled headers, and are needed to create the code file .o
# module is not dependent on any other modules, as seen here.
# SYNTAX: things to the left of ":" will be files created. Typically .o or .mod files
# For modules, both are needed, but for functions etc. only the .o file is needed.
#parameters.mod parameters.o: parameters.f95
#	$(FC) $(FCFLAGS) -c $<

%.mod %.o: %.f95
	$(FC) $(FCFLAGS) -c $<

#Ur.o: Ur.f95
#	$(FC) $(FCFLAGS) -c $<

#f.o: f.f95
#	$(FC) $(FCFLAGS) -c $<

#potential.o: potential.f95
#	$(FC) $(FCFLAGS) -c $<

# This is to make cleaning easier. Removes all the code files and header files
# It will not remove any fortran code.
.PHONY: clean
clean:
	rm -f *.o *.mod *.MOD





# Save statement for global parameters?
# Helst ikke seede med tiden? Ikke reproduserbar.
