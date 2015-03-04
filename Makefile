CC=icc
CXX=icpc

FLAGS=-g3 -parallel -openmp
OPTFLAGS=-O3 -qopt-report=3 -ip -xavx -axavx -fno-alias -fno-fnalias

CPPFLAGS= -I/home/morpheus/utils/gsl-install/include
LDFLAGS= -L/usr/lib -L/usr/local/lib -L/home/morpheus/utils/gsl-install/lib
LIBS=-lgsl -lgslcblas -lm 
LOCALINC=-I${PWD}/include/


fJmodels:
	$(CXX) -o fJmodels $(FLAGS) $(OPTFLAGS) $(CPPFLAGS) $(LOCALINC) $(LDFLAGS) fJmodels.cpp pot/*.cpp delta/*.cpp uvOrb/*.cpp df/*.cpp integ/*.cpp io/*.cpp libs/*.cpp $(LIBS)

clean:
	rm fJmodels *.optrpt config.*
