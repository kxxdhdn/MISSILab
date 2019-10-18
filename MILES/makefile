#*****************************************************************************
#*
#*       Makefile for the miles, using the static library
#*
#*****************************************************************************


# Common part
mkpath = /Users/dhu/PhD/SwING/Fortran/
mkcom = $(mkpath)make_common.mk
include $(mkcom)

# Supplementary modules
MODSS = factable.f90 \
        auxil.f90

# Programs
SRCS = $(wildcard test_*.f90)# miles.f90)
OBJSS = $(MODSS:.f90=.o)
EXEC = $(SRCS:.f90=)

# Generate all executables
all : $(EXEC)
% : $(OBJSS) %.o $(LIBDIR)$(LIBAR)
	$(FC) $(LDFLAGS) -o $@ $^ -l$(LIB) -L$(LIBDIR) -I$(LIBDIR)

# Individual file compilation
%.o : %.f90
	$(FC) $(FFLAGS) -c $? -I$(LIBDIR)
.PRECIOUS: %.o # keep intermediate object files

# House cleaning
clean :
	@echo "CLEANING THE TESTS...\n"
	@rm -f *.o *.mod
mrproper : clean
	@echo "MR. PROPERING THE TESTS...\n"
	@rm -f $(EXEC)
