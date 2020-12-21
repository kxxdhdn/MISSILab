#******************************************************************************
#*
#*            Common Data to all Makefiles of the MILES library
#*
#******************************************************************************


## Load SwING lib
##----------------
HOST = $(shell hostname)
ifeq ($(HOST),dapmcw156)
	MKSW = /Users/dhu/SwING/Fortran/make_common.mk
else ifeq ($(HOST),iclust*)
	MKSW = /dsm/herschel10/nuages/SwING/Fortran/make_common.mk
endif
$(info The current machine is $(HOST))

include $(MKSW)

## MILES modules
##---------------
MILMODS = datable.f90 \
          auxil.f90 \

## Name of the library and its paths
##-----------------------------------
MILIB = miles
MILIBAR = lib$(MILIB).a
ifeq ($(HOST),dapmcw156)
	MILIBDIR = /Users/dhu/ownCloud/MILES/
else ifeq ($(HOST),iclust*)
	MILIBDIR = /dsm/herschel1/nuages/dhu/MISSILE/MILES/
endif

## MILES Banner
##--------------
define MILBAN

>>>>>>>>>>>>>>>>>>>>>>>>>>
 MAKING THE MILES LIBRARY                     
<<<<<<<<<<<<<<<<<<<<<<<<<<

endef
$(info $(MILBAN))
