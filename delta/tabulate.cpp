/*
 * tabulate.cpp
 *
 *  Created on: Feb 22, 2015
 *      Author: L. Posti
 */

#include "Grid.h"
#include "Delta.h"
#include "Potential.h"
#include "Utils.h"

void tabulateDelta(Potential *p){
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

	for(int nE=0; nE<NGRID; nE++){//loop over E
		Egrid[nE]=Ecgrid[1]+.5*pow(Vmax[1]*(nE+1)/(double)NGRID,2);

		struct eqsPar par; par.p = *p; par.R0=Rgrid[1]; par.E=Egrid[nE]; par.Lz=.25*Lzgrid[1]; par.L=0;
		Delta * dd = new Delta (&par);
		Dgrid[nE]=dd->getBestDelta();

		delete dd;
	}

	delArr(Lzgrid); delArr(Vmax); delArr(Ecgrid); delArr(Rgrid);
}

double Deltafn(const double E){
	if(E<Egrid[0]) return sqrt(MAX(.000001,Dgrid[0]));
	else if(E<Egrid[NGRID-1]) return sqrt(MAX(.000001,linterp(Egrid,Dgrid,NGRID,E)));
	else return sqrt(MAX(.000001,Dgrid[NGRID-1]));
}
