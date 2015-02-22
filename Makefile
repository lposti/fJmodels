FLAGS_GNU   = -g3 -fopenmp
OPT_GNU     = -O3 -ftree-vectorizer-verbose=1 -march=corei7-avx -mtune=corei7-avx -mavx
FLAGS_INTEL = -g3 -parallel -openmp
OPT_INTEL   = -O3 -qopt-report=3 -ip -ipo -xavx -axavx -fno-alias -fno-fnalias
INCLUDE     = -I${PWD}/include/ -I${HOME}/utils/gsl-install/include
LIBS        = -L${HOME}/utils/gsl-install/lib
LDLIBS      = -lgsl -lgslcblas -lm

gnu:
	g++ -o fJmodels $(FLAGS_GNU) $(OPT_GNU) $(INCLUDE) $(LIBS) fJmodels.cpp pot/*.cpp delta/*.cpp uvOrb/* $(LDLIBS)


intel:
	icpc -o fJmodels $(FLAGS_INTEL) $(OPT_INTEL) $(INCLUDE) $(LIBS) fJmodels.cpp pot/*.cpp delta/*.cpp uvOrb/* $(LDLIBS)

clean:
	rm fJmodels *.optrpt
