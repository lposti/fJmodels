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

double ar[NR], Dgrid[NGRID], Egrid[NGRID];
double ** __restrict phil, ** __restrict Pr, ** __restrict Pr2,
       ** __restrict rhl, ** __restrict vrotl, ** __restrict sigRl,
       ** __restrict sigpl, ** __restrict sigzl, ** __restrict sigRzl;

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

void twoComp(struct fJParams fJP){

	double ** __restrict rhl2, ** __restrict vrotl2, ** __restrict sigRl2,
    ** __restrict sigpl2, ** __restrict sigzl2, ** __restrict sigRzl2;
	rhl2 = mat<double>(NR,NPOLY);   vrotl2 = mat<double>(NR,NPOLY);
	sigRl2 = mat<double>(NR,NPOLY); sigpl2 = mat<double>(NR,NPOLY);
	sigzl2 = mat<double>(NR,NPOLY); sigRzl2 = mat<double>(NR,NPOLY);

	SetModel(fJP,1);

	Potential p(1);
	p.selectGuessRho(fJP.modName);
	p.computeGuessRhl(); p.computePhil(p.rhlP);

	SetModel(fJP,2);
	Potential p2(2);
	p2.selectGuessRho(fJP.modName2);
	p2.computeGuessRhl(); p2.computePhil(p2.rhlP);

	/*
	 *  Now I can safely pass p since total phils and derivatives
	 *  are that of the total potential.
	 */
	tabulateDelta(&p);

	for (int k=0; k<5; k++){
		printf("\n Iter:%d\n",k);

		setDF(fJP,&p,1);
		computeNewPhi(&p,rhl,sigRl,sigpl,sigzl,sigRzl,vrotl);

		setDF(fJP,&p,2);
		computeNewPhi(&p2,rhl2,sigRl2,sigpl2,sigzl2,sigRzl2,vrotl2);
		for (int i=0; i<NR; i++) printf("%f %f %f %f %f %f %f\n",rhl[i][0],rhl2[i][0],sigRl[i][0],sigRl2[i][0],phil[i][0],p.philP[i][0],p2.philP[i][0]);
		writeOut(fJP,k,1,rhl,sigRl,sigpl,sigzl,sigRzl,vrotl,p.philP);
		writeOut(fJP,k,2,rhl2,sigRl2,sigpl2,sigzl2,sigRzl2,vrotl2,p2.philP);
	}
}

int main(int argc, char **argv){

	/*
	 * Initializing parameters
	 */

	time_t start = clock();

	phil = mat<double>(NR,NPOLY);  Pr = mat<double>(NR,NPOLY);    Pr2 = mat<double>(NR,NPOLY);
	rhl = mat<double>(NR,NPOLY);   vrotl = mat<double>(NR,NPOLY);
	sigRl = mat<double>(NR,NPOLY); sigpl = mat<double>(NR,NPOLY);
	sigzl = mat<double>(NR,NPOLY); sigRzl = mat<double>(NR,NPOLY);

	struct fJParams fJP = readParam(); printParam(fJP);

	SetGrid(50.);

	twoComp(fJP);
	if (0) {

		SetModel(fJP);

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
			//for (int i=0; i<NR; i++) printf("%f %f %f\n",ar[i],p(ar[i],0),ev_dens<double>(ar[i],0));
			computeNewPhi(&p);
			writeOut(fJP,k);
			vir2(&p);
		}
	}

	printf("\n----> Elapsed time of computation: %7.5f s\n",(clock()-start) / (double) CLOCKS_PER_SEC);
}



