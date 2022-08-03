#*****************************************************************************
#*
#*                        Makefile for the programs
#*
#*****************************************************************************


# Common part
include ../make_swing.mk

# Programs
SRCS = $(wildcard *.f90)
OBJS = $(MODS:.f90=.o)
EXEC = $(SRCS:.f90=)

# Generate all executables
all : $(EXEC) lib

# Make the library archive
lib : $(LIBAR)
$(LIBAR) : $(OBJS)
	@echo "\nCREATING THE STATIC LIBRARY...\n"
	ar cr $(LIBAR) $^

# Linking
% : $(OBJS) %.o
	$(FC) $(LDFLAGS) -o $@ $^

# Individual file compilation
%.o : ../*/%.f90
	$(FC) $(FFLAGS) -c $?
.PRECIOUS: %.o # keep intermediate object files

# House cleaning
clean :
	@echo "CLEANING THE LIBRARY...\n"
	@rm -f *.o
mrproper : clean
	@echo "MR. PROPERING THE LIBRARY...\n"
	@rm -f *.mod $(EXEC) $(LIBAR)
