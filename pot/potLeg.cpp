/*
 * potLeg.cpp
 *
 *  Created on: Feb 24, 2015
 *      Author: morpheus
 */

#include "Grid.h"
#include "Potential.h"
#include "Integ.h"
#include "Utils.h"
#include "UtilsLeg.h"
#include "Progressbar.h"

/*
 *  Compute Legendre coefficients for new rho
 *  integrating the DF
 */
void computeRhl(Potential *p){

	double **rho  = mat<double>(NR,NGAUSS);
	double **poly = mat<double>(NPOLY,NGAUSS);

	double pol[NPOLY],ci[NGAUSS],si[NGAUSS],wi[NGAUSS];
	gauleg<double>(0,1,ci,wi);

	for(int i=0; i<NGAUSS; i++){
			si[i]=sqrt(1-ci[i]*ci[i]);
			evenlegend(pol,ci[i]);
			for(int np=0; np<NPOLY; np++) poly[np][i]=2*pol[np]*wi[i];
		}

	/* initialize bar */
	ProgressBar bar(60);
	bar.init(NR*NGAUSS);

	/* integrate DF */
	for(int nr=0; nr<NR; nr++)
		for(int ng=0; ng<NGAUSS; ng++){
			bar.update(nr*NGAUSS+ng+1);
			rho[nr][ng]=rhofDF(ar[nr]*si[ng],ar[nr]*ci[ng],p);
		}

	/* finalize bar */
	bar.fillSpace("..done potential computation!!\n\n");

	/* set rhl to zeros */
	for(int nr=0; nr<NR; nr++)
		for(int np=0; np<NPOLY; np++)
			rhl[nr][np]=0.;

	/* compute rhl from new rho */
	for (int nr=0; nr<NR; nr++)
		for(int np=0; np<NPOLY; np++)
			for(int ng=0; ng<NGAUSS; ng++)
				rhl[nr][np]+=poly[np][ng]*rho[nr][ng];

	rhl[0][0]=2*rho[0][0];
	for(int i=1;i<NPOLY;i++) rhl[0][i]=0;
}

/*
 *  Compute new Potential (phil,Pr,Pr2)
 *  Computes rhl first and calls the integrator for the DF
 */
void computeNewPhi(Potential *p){

	computeRhl(p);
	if (p->canEv) p->computePhil();
}


