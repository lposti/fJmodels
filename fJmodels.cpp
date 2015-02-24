/*
 * fJmodels.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: morpheus
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Grid.h"
#include "Utils.h"
#include "readParam.h"
#include "Potential.h"
#include "models.h"
#include "Delta.h"
#include "Tabulate.h"
#include "uvOrb.h"
#include "DF.h"
#include "Integ.h"
//#include "press.h"
#include <iostream>

const double q=1,q2=q*q;
const double mass=1.,J0=1.,r0=pow(J0,2)/mass;
double ar[NR], rhl[NR][NPOLY], phil[NR][NPOLY], Pr[NR][NPOLY], Pr2[NR][NPOLY];
double Dgrid[NGRID], Egrid[NGRID];

double Hpot(double r){
	return -1./(1.+r);
}

void testInteg(Potential *p){

	FILE*fp=fopen("rho.dat","w");
	for(int n=0; n<NR; n++)
		fprintf(fp,"%f %f\n",ar[n],rhofDF(ar[n],0.,p));
}

int main(int argc, char **argv){

	/*
	 * Initializing parameters
	 */

	time_t start = clock();
	SetGrid(50.);

	struct fJParams fJP = readParam();
	std::cout << fJP.modName << " " << fJP.dphi_h_in << " " << fJP.dz_h_in << " " << fJP.dphi_g_in << " " << fJP.dz_g_in << std::endl;

	if (1) {
		Potential p;
		p.selectGuessRho(fJP.modName);
		p.computeGuessRhl(); p.computePhil();
		for (int i=0; i<NR; i++) printf("%f %f %f %f %f %f %f\n",ar[i],p.rhoGuess(ar[i],0.),rhl[i][0],phil[i][0],p(ar[i],0),p(0,ar[i]),Hpot(ar[i]));

		//struct eqsPar par; par.p = p; par.R0=1.7; par.E=-0.2; par.Lz=0.03; par.L=0;//par.E = p(ar[5],0); par.Lz=.2*ar[1]*sqrt(-2*p(ar[1],0)); par.L = par.Lz; par.R0 = ar[1];
		//Delta d(&par);
		printf("\n=============================================\n");

		printf("\n-----TEST\n");
		tabulateDelta(&p);

		setDF(fJP.dphi_h_in,fJP.dz_h_in,fJP.dphi_g_in,fJP.dz_g_in,&p);
		testInteg(&p);
	}

	printf("\n----> Elapsed time of computation: %7.5f s\n",(clock()-start) / (double) CLOCKS_PER_SEC);
}



