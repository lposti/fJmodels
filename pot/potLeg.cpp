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
void computeRhl(Potential *p, double **rhlH=rhl,double ** sigRlH=sigRl,
		double **sigplH=sigpl,double ** sigzlH=sigzl,double **sigRzlH=sigRzl,
		double **vrotlH=vrotl){

	double **rho  = mat<double>(NR,NGAUSS), **vrot  = mat<double>(NR,NGAUSS);
	double **sigR = mat<double>(NR,NGAUSS), **sigp  = mat<double>(NR,NGAUSS),
		   **sigz = mat<double>(NR,NGAUSS), **sigRz = mat<double>(NR,NGAUSS);

	double ci[NGAUSS],si[NGAUSS],wi[NGAUSS];

	gauleg<double>(0,1,ci,wi);
	for(int i=0; i<NGAUSS; i++)
		si[i]=sqrt(1-ci[i]*ci[i]);

	/* initialize bar */
	ProgressBar bar(60);
	bar.init(NR*NGAUSS);

	/* integrate DF */
	for(int nr=0; nr<NR; nr++)
		for(int ng=0; ng<NGAUSS; ng++){
			bar.update(nr*NGAUSS+ng+1);
			rho[nr][ng]=rhofDF(ar[nr]*si[ng],ar[nr]*ci[ng],p,vrot[nr]+ng,sigR[nr]+ng,sigp[nr]+ng,sigz[nr]+ng,sigRz[nr]+ng);
		}

	/* finalize bar */
	bar.fillSpace("..done potential computation!!\n\n");

	/*
	 * 03/09/14: linear interpolation in the first grid point!
	 * 			 it was causing troubles in the potential computation
	 * 			 since the density was extremely large, implying very large Phi
	 */
	rho[0][0] = (rho[2][0]-rho[1][0])/(ar[2]-ar[1])*(ar[0]-ar[1])+rho[1][0];
	get_Ylm<double>(rho,rhlH,vrot,vrotlH,sigR,sigRlH,sigp,sigplH,sigz,sigzlH,sigRz,sigRzlH);

}

/*
 *  Compute new Potential (phil,Pr,Pr2)
 *  Computes rhl first and calls the integrator for the DF
 */
void computeNewPhi(Potential *p,double **rhlH,double ** sigRlH,
		double **sigplH,double ** sigzlH,double **sigRzlH,double **vrotlH){

	computeRhl(p,rhlH,sigRlH,sigplH,sigzlH,sigRzlH,vrotlH);
	p->computePhil(rhlH);
}


