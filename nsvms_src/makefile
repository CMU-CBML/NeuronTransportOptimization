# PETSC_DIR=/ocean/projects/eng170006p/angranl/lib/petsc
PETSC_DIR=/ocean/projects/eng170006p/angranl/lib/petsc-3.15.0-gcc10.2

include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#include ${PETSC_DIR}/lib/petsc/conf/test
Compiler=mpicxx
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS}  ${PSOURCECXX}  -std=c++11 -fopenmp -g -Wall
#-D_GLIBCXX_USE_CXX11_ABI=0
all: nsvms

nsvms: BasicDataStructure.o UserSetting.o NS_3Dsteady.o UserSetting2D.o NS_2Dsteady.o Utils.o L2Projection.o main.o 
	$(Compiler) BasicDataStructure.o UserSetting.o NS_3Dsteady.o UserSetting2D.o NS_2Dsteady.o Utils.o main.o L2Projection.o -o nsvms ${PETSC_LIB} $(CFLAGS)
	
# nsvms: BasicDataStructure.o UserSetting.o NS_3Dsteady.o main.o 
	# $(Compiler) BasicDataStructure.o UserSetting.o NS_3Dsteady.o main.o -o nsvms ${PETSC_LIB} $(CFLAGS)

BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	-c BasicDataStructure.cpp $(CFLAGS)
	
UserSetting.o: UserSetting.cpp
	$(Compiler)	-c UserSetting.cpp $(CFLAGS)

UserSetting2D.o: UserSetting2D.cpp
	$(Compiler)	-c UserSetting2D.cpp $(CFLAGS)

Utils.o: Utils.cpp
	$(Compiler)	-c Utils.cpp $(CFLAGS)

L2Projection.o: L2Projection.cpp
	$(Compiler)	-c L2Projection.cpp $(CFLAGS)

NS_3Dsteady.o: NS_3Dsteady.cpp
	$(Compiler)	-c NS_3Dsteady.cpp $(CFLAGS)

NS_2Dsteady.o: NS_2Dsteady.cpp
	$(Compiler)	-c NS_2Dsteady.cpp $(CFLAGS)

main.o: main.cpp
	$(Compiler)	-c main.cpp $(CFLAGS)

clean:
	rm *o nsvms