# Simple Makefile for a Fortran program
# Program: poisson.f90
# Module: precision.f90 (used by poisson.f90)

# Compiler and flags
# Fortran compiler
FC = gfortran
# Optimization level 3
FFLAGS = -O3

# Define source files and the final executable name
MODULE = precision.f90
MAIN = poisson.f90
EXEC = poisson

# Define object files: these are compiled versions of the source files
OBJS = precision.o poisson.o

# Default target: builds the executable
# This rule says: to build `poisson`, first make sure all object files are compiled
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS)

# Rule to compile the module
# (modules must be compiled before the files that use them)
precision.o: precision.f90
	$(FC) $(FFLAGS) -c precision.f90

# Rule to compile the main program
# Depends on both poisson.f90 and precision.o
poisson.o: poisson.f90 precision.o
	$(FC) $(FFLAGS) -c poisson.f90

# Utility target: clean up compilation artifacts
# Run `make clean` to remove object files, module files, and the executable
clean:
	rm -f *.o *.mod $(EXEC)
