/*
 * print.cpp
 *
 *  Created on: Oct 22, 2014
 *      Author: LP
 */

#include <stdio.h>
#include <stdlib.h>
#include "print.h"



#ifdef PRINTDFH

void openDFH(){ dfhp = fopen("DFen.dat","w"); }

void printDFH(double Jr, double Jphi, double Jz, double DF, double E){
	fprintf(dfhp, "%e %e %e %e %e\n",Jr,Jphi,Jz,DF,E);
}

void closeDFH(){ fclose(dfhp); }

#endif
