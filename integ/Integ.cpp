/*
 * integ.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: morpheus
 */

#include <math.h>
#include "DF.h"
#include "Potential.h"
#include "GaussQuad.h"

//double rhoInteg(const double v[], const double x[], const double Ve){

struct IntegPar { double * x; double Ve;};
double rhoInteg(const double vr, const double vphi, const double vz, void * params){

	struct IntegPar * par = (struct IntegPar*) params;
	double Ve = par->Ve;
	double *x = par->x;

	/* rescaling the integral to [0,1]x[-1,1]x[0,1] */
	double V[3] = {1e-6+Ve*vr, 1e-6+Ve*vphi, 1e-6+Ve*vz};
	double jacob = Ve*Ve*Ve;

	return jacob*df(x,V);
}

/*
 *  rho from DF integration
 */
double rhofDF(double R, double z, Potential *p){
	double Phi_h=(*p)(R,z);
	double Rz[3]={R,z,Phi_h},Ve=sqrt(-2*(Phi_h-(*p)(100,100)));

	struct IntegPar par = {&Rz[0],Ve};
	return 4*Int3D_011101(&rhoInteg,&par);
}


