 
CXX=g++-9
LIBS= -llapack -lblas   -lgfortran 
INCLUDE= -I/usr/local/Cellar/eigen/3.3.7/include/eigen3/Eigen/src/misc/


all:
	gfortran -g1 -c utilmod2.f90 diagev.f90 
	gfortran -g1 -c commutp.f90 
	gfortran -g1 -c timing.f90 dochol.f90 
	gfortran -g1 -c initmd.f90  eldyn.f90 

	${CXX}    -g1 -c main.cpp                  ${INCLUDE}
	${CXX}    -g1 -c diagev2.cpp               ${INCLUDE}
	${CXX}    -g1 -c commutp2.cpp              ${INCLUDE}
	${CXX}    -g1 -c eldyn_ort.cpp             ${INCLUDE}
	${CXX}    -g1 -c dochol2.cpp               ${INCLUDE}
	${CXX}    -g1 -c eldyn_nonort.cpp          ${INCLUDE}
	${CXX}    -g1 -c test-main-propagator.cpp  ${INCLUDE}
	#-------
	${CXX}  -g1 -o test-propagator test-main-propagator.o eldyn_nonort.o eldyn.o commutp.o commutp2.o  diagev.o  diagev2.o  utilmod2.o dochol2.o dochol.o timing.o  ${INCLUDE}  ${LIBS} 

        


