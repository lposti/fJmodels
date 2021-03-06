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
	FILE * ff=fopen("pottest.dat","w");
	FILE * fp=fopen("potMN.dat","w");
	FILE * fd=fopen("dtest.dat","w");
	double * gg = (double *) malloc(100*sizeof(double));
	gg=giveLinGrid<double>(0.01,10,100);

	for (int i=0; i<100; i++){
		for (int j=0; j<100; j++)
			fprintf(ff,"%e ", p(gg[i],gg[j]));
		fprintf(ff,"\n");
	}
	for (int i=0; i<100; i++){
		for (int j=0; j<100; j++)
			fprintf(fp,"%e ", -1./sqrt(gg[i]*gg[i] + pow(1.+sqrt(gg[j]*gg[j]+0.3*0.3),2))  );
		fprintf(fp,"\n");
	}
	for (int i=0; i<100; i++){
		for (int j=0; j<100; j++)
			fprintf(fd,"%e ", ev_dens(gg[i],gg[j]));
		fprintf(fd,"\n");
	}
	fclose(ff); fclose(fd); fclose(fp);
	exit(1);
	*/
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

	SetGrid(200.);
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



