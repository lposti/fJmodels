/*
 * fJmodels.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: morpheus
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Grid.h"
#include "Utils.h"
#include "Potential.h"
#include "models.h"
#include "Delta.h"
#include "Tabulate.h"
#include "uvOrb.h"
//#include "press.h"

const double q=1,q2=q*q;
const double mass=1.,J0=1.,r0=pow(J0,2)/mass;
double ar[NR], rhl[NR][NPOLY], phil[NR][NPOLY], Pr[NR][NPOLY], Pr2[NR][NPOLY];
double Dgrid[NGRID], Egrid[NGRID];

double Hpot(double r){
	return -1./(1.+r);
}

int main(int argc, char **argv){

	/*
	 * Initializing parameters
	 */

	time_t start = clock();
	SetGrid(50.);
	if (1) {
		Potential p;
		p.selectGuessRho("Hernquist");
		p.computeGuessRhl(); p.computePhil();
		for (int i=0; i<NR; i++) printf("%f %f %f %f %f %f %f\n",ar[i],p.rhoGuess(ar[i],0.),rhl[i][0],phil[i][0],p(ar[i],0),p(0,ar[i]),Hpot(ar[i]));

		struct eqsPar par; par.p = p; par.R0=1.7; par.E=-0.2; par.Lz=0.03; par.L=0;//par.E = p(ar[5],0); par.Lz=.2*ar[1]*sqrt(-2*p(ar[1],0)); par.L = par.Lz; par.R0 = ar[1];
		Delta d(&par);
		printf("\n=============================================\n");

		//printf("%f %f\n",d.getBestDelta(),get_closed(d.R0,par.Lz,par.E,&p));

		printf("\n-----TEST\n");
		tabulateDelta(&p);

		double *x=arr<double>(3),*v=arr<double>(3);
		x[0]=ar[NR-3]; x[1]=ar[0]; x[2]=p(x[0],x[1]); v[0]=sqrt(-2*(x[2]-p(100,100)))/5; v[1]=sqrt(-2*(x[2]-p(100,100)))/10; v[2]=sqrt(-2*(x[2]-p(100,100)))/20;

		double R=x[0],Lz=R*v[1],Phigl=x[2];
		double H=.5*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2))+Phigl,DD=Deltafn(H);
		uvOrb uv(DD,Lz,Phigl,x,v,&p);
		printf("Ju=%e Jv=%e\n",uv.Ju(),uv.Jv());
	}

	printf("\n----> Elapsed time of computation: %7.5f s\n",(clock()-start) / (double) CLOCKS_PER_SEC);
}



