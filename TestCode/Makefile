#PETSC_DIR=/mnt/c/Users/LAR/Documents/GroupCode/lib/petsc
# PETSC_DIR=~/Documents/petsc
PETSC_DIR=/pylon5/eg560mp/angranl/lib/petsc

# CLFAGS          =
# FFLAGS	        =
# CPPFLAGS        =
# FPPFLAGS        =
# LOCDIR          = src/ts/examples/tutorials/advection-diffusion-reaction/
# EXAMPLESC       = ex1.c ex2.c ex3.c ex4.c ex5.c ex5opt_ic.c
# EXAMPLESF       = shashi.F90
# EXAMPLESFH      =
# MANSEC          = TS
# DIRS            =
# CLEANFILES      =  SA-data/*

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

# Compiler=mpiicpc
# Compiler=pgcc
Compiler=mpiicc
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS}  ${PSOURCECXX} -Wall
#-D_GLIBCXX_USE_CXX11_ABI=0
all: main

main: ex5opt_ic.o
	$(Compiler) ex5opt_ic.o -o main ${PETSC_LIB} $(CFLAGS)

ex5opt_ic.o: ex5opt_ic.c
	$(Compiler)	-c ex5opt_ic.c $(CFLAGS)

# clean:
# 	rm *o