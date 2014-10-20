CFL	= -g -fopenmp -O3 -funsafe-math-optimizations  -msse4.1 -Wstrict-aliasing=2 -march=corei7 -mtune=corei7

all: fJmodels

.cpp.o:
	g++ -c $(CFL) $*.cpp

pmongo_c.a:	pmongo_c.cpp
	g++ -c $(CFL) $(OPT) pmongo_c.cpp
	test -e pmongo_c.a  && rm pmongo_c.a || echo not there
	ar -qs pmongo_c.a pmongo_c.o

press_c.a: pressstuff.cpp
	g++ -c $(CFL) $(OPT) pressstuff.cpp #-I/usr/include/x86_64-linux-gnu/c++/4.8
	test -e press_c.a && rm press_c.a || echo not there
	ar -qs press_c.a pressstuff.o
# to make it with g++ ln -s from /usr/include/x86_64-linux-gnu/c++/4.8 to /usr/include/c++/4.8 the required files

leg_pot2.o: leg_pot2.cpp

potLEG.o: potLEG.cpp leg_pot2.h

isopot.o: isopot.cpp

satoh.o: satoh.cpp

uv_orb3.o: uv_orb3.cpp uv_orb3.h Rzuv.o leg_pot2.h

Rzuv.o:	Rzuv.cpp

moments.o: moments.cpp oct_int.h

tables3.o: tables3.cpp grid3.h uv_orb.h

node1.o: node1.cpp node1.h

node2.o: node2.cpp node2.h

node3.o: node3.cpp node3.h

oct_int.o: oct_int.cpp node1.h node2.h node3.h

isodf4.o: isodf4.cpp

plot_iso: plot_iso.cpp potLEG.o leg_pot2.o isopot.o pmongo_c.a press_c.a
	g++ $(CFL) plot_iso.cpp leg_pot2.o potLEG.o isopot.o pmongo_c.a press_c.a -o plot_iso

fJmodels: fJmodels.cpp get_closed.o getDelta.o potLEG.o leg_pot2.o isopot.o Rzuv.o uv_orb3.o tables3.o node1.o node2.o node3.o oct_int.o moments.o isodf4.o check_acts.o pmongo_c.a press_c.a
	g++ $(CFL) test3.cpp get_closed.o getDelta.o potLEG.o leg_pot2.o isopot.o Rzuv.o uv_orb3.o tables3.o node1.o node2.o node3.o oct_int.o moments.o isodf4.o check_acts.o pmongo_c.a press_c.a -o fJmodels

test_tab: test_tab.cpp get_closed.o getDelta.o potLEG.o leg_pot2.o isopot.o Rzuv.o uv_orb2.o tables3.o node1.o node2.o node3.o oct_int.o moments.o isoDF2.o pmongo_c.a press_c.a
	g++ $(CFL) test_tab.cpp get_closed.o getDelta.o potLEG.o leg_pot2.o isopot.o Rzuv.o uv_orb2.o press_c.a tables3.o node1.o node2.o node3.o oct_int.o moments.o isoDF2.o pmongo_c.a -o test_tab

testpot: testpot.cpp potLEG.o leg_pot2.o isopot.o
	g++ $(CFL) testpot.cpp potLEG.o leg_pot2.o isopot.o -lpmongo_c -lpress_c -o testpot

rowley: rowley.cpp potLEG.o leg_pot2.o potLEG.h leg_pot2.h
	g++ $(CFL) rowley.cpp potLEG.o leg_pot2.o -lpmongo_c -lpress_c -o rowley

clean:
	rm -rf *.o relax plot_iso

veryclean:
	rm -rf *.o *.a relax plot_iso
