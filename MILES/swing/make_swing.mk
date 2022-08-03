#*****************************************************************************
#*
#*                        Makefile for the programs
#*
#*****************************************************************************

## Compiler
##----------
FC = h5fc
FFLAGS = -Wall -O3 -g # Standard (mono-processor)
# FFLAGS = -Wall -O3 -g -fcheck=all # debug
#FFLAGS = -Wall -O3 -ftree-parallelize-loops=4 # Mac (parallelization)
#FFLAGS = -warn all -O3 -heap-arrays -traceback # irfucoast intel compiler
#FFLAGS = -warn all -check all -g -debug -heap-arrays -traceback # debug
LDFLAGS = $(FFLAGS)

## SwING modules
##---------------
MODS = utilities.f90 \
       constants.f90 \
       arrays.f90 \
       inout.f90 \
       fft_specials.f90 \
       linear_system.f90 \
       matrices.f90 \
       interpolation.f90 \
       integration.f90 \
       special_functions.f90 \
       statistics.f90 \
       distributions.f90 \
       chi2_minimization.f90 \
       adaptative_grid.f90 \
       random.f90 \
       statistical_physics.f90 \
       grain_optics.f90 \

# Name of the library and its paths
##-----------------------------------
LIB = swing
LIBAR = lib$(LIB).a
HOST = $(shell hostname)
ifeq ($(HOST),dapmcw18)
       LIBDIR = /Users/dhu/Github/MISSILE/MILES/swing/Programs/
else
       LIBDIR = /dsm/herschel1/nuages/Codes/MILES/swing/Programs/
endif

## SwING banner
##--------------
define BANNER

*******************************************************************************
                      MAKING SwING LIBRARY snippet
*******************************************************************************

endef
$(info $(BANNER))
