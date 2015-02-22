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
//#include "press.h"

const double q=1,q2=q*q;
const double mass=1.,J0=1.,r0=pow(J0,2)/mass;
double ar[NR], rhl[NR][NPOLY], phil[NR][NPOLY], Pr[NR][NPOLY], Pr2[NR][NPOLY];

double Hpot(double r){
	return -1./(1.+r);
}

void testDelta(Potential *p){
	double * Lzgrid = arr<double>(NR), * Vmax = arr<double>(NR), * Ecgrid = arr<double>(NR), * Rgrid = arr<double>(NR);
	double Rmax=ar[NR-1];
	double Lzmax=Rmax*0.98*sqrt(-2*(*p)(Rmax,0)), Rmin=1e-3*Rmax;
	double Lzmin=Rmin*sqrt(Rmin*p->dR(Rmin,0)),R=.01*Rmax;

	for(int i=0; i<NR; i++){
		Lzgrid[i]=Lzmin+i*(Lzmax-Lzmin)/(double)(NR-1);
		double Lzsq=pow(Lzgrid[i],2); R=Rgrid[i]=GetRc<double>(Lzsq,R,p);
		double vc2=R*p->dR(R,0);
		Ecgrid[i]=.5*vc2+(*p)(R,0);//Ecirc at R
		Vmax[i]=.98*sqrt(-2*Ecgrid[i]);//Range of Egrid at R

	}

	for(int nE=0; nE<NGRID*(NR-1)/NR; nE++){//loop over E
		double E=Ecgrid[1]+.5*pow(Vmax[1]*(nE+1)/(double)(NGRID*(NR-1)/NR),2);

		struct eqsPar par; par.p = *p; par.R0=Rgrid[1]; par.E=E; par.Lz=.25*Lzgrid[1]; par.L=0;
		Delta * dd = new Delta (&par);
		//printf("(E,LZ)=(%f,%f) Mine:%f Binn:%f\n",E,par.Lz,dd->getBestDelta(),get_closed(RcE(E,Rgrid[1],p),.25*Lzgrid[1],E,p));
		printf("(E,LZ)=(%f,%f) Mine:%f Rc:%f\n",E,par.Lz,dd->getBestDelta(),RcE<double>(E,1,p));

		delete dd;
	}
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
		testDelta(&p);
	}

	printf("\n----> Elapsed time of computation: %7.5f s\n",(clock()-start) / (double) CLOCKS_PER_SEC);
}



