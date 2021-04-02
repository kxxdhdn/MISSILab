#******************************************************************************
#*
#*                     Makefile for the programs
#*
#******************************************************************************


## Common part
##-------------
include ../make_miles.mk

## Tests
##-------
SRCS = $(wildcard *.f90)
OBJS = $(MILMODS:.f90=.o)
EXEC = $(SRCS:.f90=)

## Generate all executables
##--------------------------
all : $(EXEC)

## Linking
# ##-------
% : $(OBJS) %.o $(LIBDIR)$(LIBAR)
	$(FC) $(LDFLAGS) -o $@ $^ -l$(LIB) -L$(LIBDIR) -I$(LIBDIR)

## Individual file compilation
##-----------------------------
%.o : ../*/%.f90
	$(FC) $(FFLAGS) -c $? -I$(LIBDIR)
.PRECIOUS: %.o # keep intermediate object files

## House cleaning
##----------------
clean :
	@echo "CLEANING...\n"
	@rm -f *.o *.mod
mrproper : clean
	@echo "MR. PROPERING...\n"
	@rm -f $(EXEC)
