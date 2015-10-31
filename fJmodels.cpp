/*
 * fJmodels.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: L. Posti
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
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

FILE * __restrict f_lz;
/*
 * Compute 2-component models
 */
void twoComp(struct fJParams fJP){

	double ** __restrict rhl2, ** __restrict vrotl2, ** __restrict sigRl2,
    ** __restrict sigpl2, ** __restrict sigzl2, ** __restrict sigRzl2;
	rhl2 = mat<double>(NR,NPOLY);   vrotl2 = mat<double>(NR,NPOLY);
	sigRl2 = mat<double>(NR,NPOLY); sigpl2 = mat<double>(NR,NPOLY);
	sigzl2 = mat<double>(NR,NPOLY); sigRzl2 = mat<double>(NR,NPOLY);

	printf("\nSetting up two-component f(J) model...\n");
	SetModel(fJP,1);

	Potential p(1);
	p.selectGuessRho( fJP.modName!="null" ? fJP.modName : "Hernquist" );
	p.computeGuessRhl(); p.computePhil(p.rhlP);

	SetModel(fJP,2);
	Potential p2(2);
	p2.selectGuessRho( fJP.modName2!="null" ? fJP.modName2 : "Hernquist" );
	p2.computeGuessRhl(); p2.computePhil(p2.rhlP);

	updatePhil(&p, &p2);		// update the total potential
	/*
	 *  Now I can safely pass p since total phils and derivatives
	 *  are that of the total potential.
	 */
	tabulateDelta(&p);

	for (int i=0; i<NR; i++) printf("%f %f %f %f %f \n",p.rhlP[i][0],p2.rhlP[i][0],phil[i][0],p.philP[i][0],p2.philP[i][0]);
	for (int k=0; k<fJP.itermax; k++){
		printf("\n Iter:%d\n",k);

		setDF(&fJP,&p,1);
		computeNewPhi(&p,rhl,sigRl,sigpl,sigzl,sigRzl,vrotl);

		setDF(&fJP,&p,2);
		computeNewPhi(&p2,rhl2,sigRl2,sigpl2,sigzl2,sigRzl2,vrotl2);

		updatePhil(&p, &p2);		// update the total potential

		//for (int i=0; i<NR; i++) printf("%f %f %f %f %f %f %f\n",rhl[i][0],rhl2[i][0],sigRl[i][0],sigRl2[i][0],phil[i][0],p.philP[i][0],p2.philP[i][0]);
		writeOut(fJP,k,1,rhl,sigRl,sigpl,sigzl,sigRzl,vrotl,p.philP,p.PrP,p.Pr2P);
		writeOut(fJP,k,2,rhl2,sigRl2,sigpl2,sigzl2,sigRzl2,vrotl2,p2.philP,p2.PrP,p2.Pr2P);
	}
}

/*
 * Compute 1-component models
 */
void oneComp(struct fJParams fJP){

	printf("\nSetting up one-component f(J) model...\n");
	SetModel(fJP);

	Potential p;
	p.selectGuessRho( fJP.modName!="null" ? fJP.modName : "Hernquist" );
	p.computeGuessRhl(); p.computePhil();

/*
	Potential ext(2);
	ext.selectGuessRho("NFWext");
	ext.computeGuessRhl(); ext.computePhil();
	updatePhil(&p, &ext);		// update the total potential
*/

	tabulateDelta(&p);

	setDF(&fJP,&p,1);

	for (int k=0; k<fJP.itermax; k++){
		printf("\n Iter:%d\n",k);
		f_lz = fopen("f_lz.dat","w");

		//double **philOLD=mat<double>(NR,NPOLY), **PrOLD=mat<double>(NR,NPOLY),**Pr2OLD=mat<double>(NR,NPOLY);
		//savePhi(philOLD,PrOLD,Pr2OLD);
		computeNewPhi(&p,rhl,sigRl,sigpl,sigzl,sigRzl,vrotl);
		//updatePhil(&p, &ext);		// update the total potential
		//mergePhi(philOLD,PrOLD,Pr2OLD,.25);
		writeOut(fJP,k,&p);
		vir2(&p);
		fclose(f_lz);
	}

	/*
	double ci[6],si[6],wi[6];
	// write f_lz file
	gauleg<double>(0,1,ci,wi);
	for(int i=0; i<6; i++)
		si[i]=sqrt(1-ci[i]*ci[i]);

	f_lz = fopen("f_lz.dat","w");
	double Ve;
	int nn=20;
	double r,vr,vphi,vz,V[3],X[3],Phi_h;
	for(int nr=0; nr<NR; nr++){
		printf("%d\n",nr);
		r = ar[0] + 15. * float(nr)/float(NR-1);
		for(int ng=0; ng<6; ng++){
			Phi_h = p(r*si[ng],r*ci[ng]);
			X[0]=r*si[ng]; X[1]=r*ci[ng]; X[2]=Phi_h;

#pragma omp parallel
{
#pragma omp parallel for private(Ve,vr,vphi,vz)
			for (int ir=0; ir<nn; ir++){
				Ve = sqrt(-2*(Phi_h-p(10 * ar[NR-1],10 * ar[NR-1])));
				vr = 1e-6+Ve * float(ir)/float(nn-1) * .5;
				for (int ip=0; ip<nn; ip++){
					vphi =(1e-6-Ve + 2. * Ve * float(ip)/float(nn-1)) * .5;
					for (int iz=0; iz<nn; iz++){
						vz = 1e-6+Ve * float(iz)/float(nn-1) * .5;
						V[0]=vr; V[1]=vphi; V[2]=vz;
						df(X,V);
					}
				}
			}
}
		}
	}
	*/
}

int main(int argc, char **argv){

	/*
	 * Initializing parameters
	 */

	//time_t start = clock();
	double   start = omp_get_wtime();

	phil = mat<double>(NR,NPOLY);  Pr = mat<double>(NR,NPOLY);    Pr2 = mat<double>(NR,NPOLY);
	rhl = mat<double>(NR,NPOLY);   vrotl = mat<double>(NR,NPOLY);
	sigRl = mat<double>(NR,NPOLY); sigpl = mat<double>(NR,NPOLY);
	sigzl = mat<double>(NR,NPOLY); sigRzl = mat<double>(NR,NPOLY);

	/*
	 *  read parameter file:
	 *  if the program was launched with an argument interpret it as an input file.
	 *  if nothig was specified by default the program will search for a file called "param.txt"
	 *  in the root dir
	 */
	struct fJParams fJP;
	if (argc>1) fJP = readParam(argv[1]);
	else fJP = readParam();
	printParam(&fJP);

	SetGrid(100.);
	printf("=");
	printf("\n======================================================================\n");

	if (components == 1) oneComp(fJP);
	else if (components == 2) twoComp(fJP);
	else
	{
		printf("ERROR: Cannot initialize neither one- nor two-component models!");
		exit(6);
	}

	// printf("\n----> Elapsed time of computation: %7.5f s\n",(clock()-start) / (double) CLOCKS_PER_SEC);
	printf("\n----> Elapsed time of computation: %7.5f s\n",(omp_get_wtime()-start));
}



