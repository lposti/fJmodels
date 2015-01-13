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
static FILE* dfhp;


void openDFH(){ dfhp = fopen("DFen.dat","w"); }

void printDFH(double Jr, double Jphi, double Jz, double DF, double E){
	fprintf(dfhp, "%e %e %e %e %e\n",Jr,Jphi,Jz,DF,E);
}

void closeDFH(){ fclose(dfhp); }

#endif

#ifdef PRINTXVJ
static FILE* xvjp;

void openXVJ(){ xvjp = fopen("xvJ.dat","w"); }

void printXVJ(double *x, double *v, double Jr, double Jphi, double Jz){
	fprintf(xvjp, "%e %e %e %e %e %e %e %e\n",x[0],x[1],v[0],v[1],v[2],Jr,Jphi,Jz);
}

void closeXVJ(){ fclose(xvjp); }

#endif
