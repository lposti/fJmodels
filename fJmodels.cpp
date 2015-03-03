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
#include "UtilsLeg.h"
#include "readParam.h"
#include "Potential.h"
#include "models.h"
#include "Delta.h"
#include "Tabulate.h"
#include "uvOrb.h"
#include "DF.h"
#include "Integ.h"
#include "potLeg.h"
#include "GaussPts.h"
#include "Stats.h"
#include "writeOut.h"

double ar[NR], phil[NR][NPOLY], Pr[NR][NPOLY], Pr2[NR][NPOLY];
double Dgrid[NGRID], Egrid[NGRID];
double ** __restrict rhl,** __restrict vrotl,** __restrict sigRl,
       ** __restrict sigpl,** __restrict sigzl,** __restrict sigRzl;

double Hpot(double r){
	return -1./(1.+r);
}

void testInteg(Potential *p){

	double vrot,sigR,sigp,sigz,sigRz;
	FILE*fp=fopen("rho.dat","w");
	for(int n=0; n<NR; n++){
		double rhoh=rhofDF(ar[n],0.,p,&vrot,&sigR,&sigp,&sigz,&sigRz);
		fprintf(fp,"%f %f %f %f %f %f %f\n",ar[n],rhoh,vrot,sigR,sigp,sigz,sigRz);
	}
}

int main(int argc, char **argv){

	/*
	 * Initializing parameters
	 */

	time_t start = clock();

	rhl = mat<double>(NR,NPOLY);   vrotl = mat<double>(NR,NPOLY);
	sigRl = mat<double>(NR,NPOLY); sigpl = mat<double>(NR,NPOLY);
	sigzl = mat<double>(NR,NPOLY); sigRzl = mat<double>(NR,NPOLY);

	struct fJParams fJP = readParam(); printParam(fJP);

	SetGrid(50.); SetModel(fJP);

	if (1) {
		Potential p;
		p.selectGuessRho(fJP.modName);
		p.computeGuessRhl(); p.computePhil();
		//for (int i=0; i<NR; i++) printf("%f %f %f\n",ar[i],p(ar[i],0),ev_dens<double>(ar[i],0));

		printf("\n======================================================================\n");

		printf("\n-----TEST %d\n",QUADORD);
		tabulateDelta(&p);

		setDF(fJP,&p);
		//testInteg(&p);

		for (int k=0; k<5; k++){
			printf("\n Iter:%d\n",k);
			//Potential ext(false);
			//ext.selectGuessRho("NFWext");
			//ext.computeGuessRhl(); ext.computePhil();
			for (int i=0; i<NR; i++) printf("%f %f %f\n",ar[i],p(ar[i],0),ev_dens<double>(ar[i],0));
			computeNewPhi(&p);
			writeOut(fJP,k);
			vir2(&p);
		}
	}

	printf("\n----> Elapsed time of computation: %7.5f s\n",(clock()-start) / (double) CLOCKS_PER_SEC);
}



