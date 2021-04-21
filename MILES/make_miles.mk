#******************************************************************************
#*
#*            Common Data to all Makefiles of the MILES library
#*
#******************************************************************************


## Load SwING lib
##----------------
HOST = $(shell hostname)
ifeq ($(HOST),dapmcw18)
	MKSW = /Users/dhu/SwING/Fortran/make_common.mk
else
	MKSW = /dsm/herschel1/nuages/Codes/SwING/Fortran/make_common.mk
endif
$(info The current machine is $(HOST))

include $(MKSW)

## MILES modules
##---------------
MILMODS = datable.f90 \
          auxil.f90 \
          chi2_kit.f90 \
          HB_kit.f90 \

## Name of the library and its paths
##-----------------------------------
MILIB = miles
MILIBAR = lib$(MILIB).a
ifeq ($(HOST),dapmcw18)
	MILIBDIR = /Users/dhu/ownCloud/MILES/
else
	MILIBDIR = /dsm/herschel1/nuages/dhu/MISSILE/MILES/
endif

## MILES Banner
##--------------
define MILBAN

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                         MAKING THE MILES LIBRARY                     
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

endef
$(info $(MILBAN))
