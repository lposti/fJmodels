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
#include "potLeg.h"
#include <iostream>

const double q=1,q2=q*q;
const double mass=1.,J0=1.,r0=pow(J0,2)/mass;
double ar[NR], rhl[NR][NPOLY], phil[NR][NPOLY], Pr[NR][NPOLY], Pr2[NR][NPOLY];
double Dgrid[NGRID], Egrid[NGRID];

double Hpot(double r){
	return -1./(1.+r);
}

void testInteg(Potential *p){

	double sigR,sigp,sigz,sigRz;
	FILE*fp=fopen("rho.dat","w");
	for(int n=0; n<NR; n++){
		double rhoh=rhofDF(ar[n],0.,p,&sigR,&sigp,&sigz,&sigRz);
		fprintf(fp,"%f %f %f %f %f %f\n",ar[n],rhoh,sigR,sigp,sigz,sigRz);
	}
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
		for (int i=0; i<NR; i++) printf("%f %f\n",ar[i],phil[i][0]);

		printf("\n=============================================\n");

		printf("\n-----TEST\n");
		tabulateDelta(&p);

		setDF(fJP.dphi_h_in,fJP.dz_h_in,fJP.dphi_g_in,fJP.dz_g_in,&p);
		testInteg(&p);
		//computeNewPhi(&p);
		//for (int i=0; i<NR; i++) printf("%f %f\n",ar[i],phil[i][0]);
	}

	printf("\n----> Elapsed time of computation: %7.5f s\n",(clock()-start) / (double) CLOCKS_PER_SEC);
}



