#******************************************************************************
#*
#*            Common Data to all Makefiles of the MILES library
#*
#******************************************************************************


## Load SwING lib
##----------------
HOST = $(shell hostname)
ifeq ($(HOST),dapmcw18)
# 	MKSW = /Users/dhu/Github/MISSILE/MILES/swing/make_swing.mk # extracted SwING
	MKSW = /Users/dhu/SwING/Fortran/make_common.mk
else
# 	MKSW = /dsm/herschel1/nuages/dhu/MILES/swing/make_swing.mk # extracted SwING
	MKSW = /dsm/herschel1/nuages/Codes/SwING/Fortran/make_common.mk
endif
$(info The current machine is $(HOST))

include $(MKSW)

## MILES modules
##---------------
MILMODS = auxil.f90 \
          core.f90 \
          ext_chi2.f90 \
          ext_hb.f90 \

## Name of the library and its paths
##-----------------------------------
MILIB = miles
MILIBAR = lib$(MILIB).a
ifeq ($(HOST),dapmcw18)
	MILIBDIR = /Users/dhu/Github/MISSILE/MILES/
else
	MILIBDIR = /dsm/herschel1/nuages/dhu/MILES/
endif

## MILES banner
##--------------
define MILBAN

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                         MAKING THE MILES LIBRARY                     
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

endef
$(info $(MILBAN))
