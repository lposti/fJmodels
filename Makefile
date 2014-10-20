FLAGS = -g3 -fopenmp
#OPT   = -O3 -funsafe-math-optimizations  -msse4.1 -Wstrict-aliasing=2 -march=native -mtune=native
OPT    = -O2 -march=native -mtune=native
INCLUDE = -I${HOME}/utils/gsl-install/include
LIBS    = -L${HOME}/utils/gsl-install/lib
LDLIBS  = -lgsl -lgslcblas -lm

fJmodels:
	g++ -o fJmodels $(FLAGS) $(OPT) $(INCLUDE) $(LIBS) *.cpp $(LDLIBS)

clean:
	rm fJmodels
