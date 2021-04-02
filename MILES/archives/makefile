#*****************************************************************************
#*
#*             Makefile for miles, using the static library
#*
#*****************************************************************************


## Common part
##-------------
include ./make_miles.mk

## Programs
##----------
SRCS = $(wildcard *.f90)
OBJS = $(MILMODS:.f90=.o)
EXEC = $(SRCS:.f90=)

## Generate all executables
##--------------------------
all : $(EXEC) lib

## Make the library archive
##--------------------------
lib : $(MILIBAR)
$(LIBAR) : $(OBJS)
	@echo "\nCREATING THE STATIC LIBRARY...\n"
	ar cr $(MILIBAR) $^

## Linking
# ##-------
% : $(OBJS) %.o #$(LIBDIR)$(LIBAR)
	$(FC) $(LDFLAGS) -o $@ $^ -l$(LIB) -L$(LIBDIR) -I$(LIBDIR)

## Individual file compilation
##-----------------------------
%.o : %.f90
	$(FC) $(FFLAGS) -c $? -I$(LIBDIR)
.PRECIOUS: %.o # keep intermediate object files

## House cleaning
##----------------
clean :
	@echo "CLEANING THE TESTS...\n"
	@rm -f *.o *.mod
mrproper : clean
	@echo "MR. PROPERING THE TESTS...\n"
	@rm -f $(EXEC)
