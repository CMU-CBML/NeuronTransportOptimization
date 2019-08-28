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

Compiler=mpiicpc
# Compiler=pgcc
# Compiler=mpiicc
# CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS}  ${PSOURCECXX} -Wall
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS}  ${PSOURCECXX}  -std=c++11 -qopenmp -Wall
#-D_GLIBCXX_USE_CXX11_ABI=0

all: neuron_opt

neuron_opt: BasicDataStructure.o UserSetting2D.o TransportOpt2D.o main.o 
	$(Compiler) BasicDataStructure.o UserSetting2D.o TransportOpt2D.o main.o -o neuron_opt ${PETSC_LIB} $(CFLAGS)
	
BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	-c BasicDataStructure.cpp $(CFLAGS)
	
UserSetting2D.o: UserSetting2D.cpp
	$(Compiler)	-c UserSetting2D.cpp $(CFLAGS)

TransportOpt2D.o: TransportOpt2D.cpp
	$(Compiler)	-c TransportOpt2D.cpp $(CFLAGS)

main.o: main.cpp
	$(Compiler)	-c main.cpp $(CFLAGS)

# clean:
# 	rm *o neuron_opt
