/*
 * integ.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: L. Posti
 */

#include <math.h>
#include "DF.h"
#include "Potential.h"
#include "GaussQuad.h"
#include "Utils.h"


struct IntegPar { double * x; double Ve;};
double rhoInteg (const double vr, const double vphi, const double vz, void * params){

	struct IntegPar * par = (struct IntegPar*) params;
	double Ve = par->Ve;
	double *x = par->x;

	/* rescaling the integral to [0,1]x[-1,1]x[0,1] */
	double V[3] = {1e-6+Ve*vr, 1e-6+Ve*vphi, 1e-6+Ve*vz};
	double jacob = Ve*Ve*Ve;

	return jacob*df(x,V);
}

/* Version with also sigma integration: returns a vector of dimension 5 */
struct vec6d rhoInteg_vec (const double vr, const double vphi, const double vz, void * params){

	struct IntegPar * par = (struct IntegPar*) params;
	double Ve = par->Ve;
	double *x = par->x;

	/* rescaling the integral to [0,1]x[-1,1]x[0,1] */
	double V[3] = {1e-6+Ve*vr, 1e-6+Ve*vphi, 1e-6+Ve*vz};
	double jacob = Ve*Ve*Ve;
	double DF = df(x,V);

	struct vec6d res = {jacob*DF,jacob*DF*V[1],jacob*DF*V[0]*V[0],
						jacob*DF*V[1]*V[1],jacob*DF*V[2]*V[2],
						jacob*DF*V[0]*V[2] };

	return res;
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

/*
 *  rho from DF integration: with sigma(R,z,phi,Rz) computed
 */
double rhofDF(double R, double z, Potential *p, double * vrot, double * sigR,
		      double * sigp, double * sigz, double * sigRz){
	double Phi_h=(*p)(R,z);
	double Rz[3]={R,z,Phi_h},Ve=sqrt(-2*(Phi_h-(*p)(100,100)));

	struct IntegPar par = {&Rz[0],Ve};
	double * out = arr<double>(6);
	double dens=4*Int3D_011101_vec(&rhoInteg_vec,&par,out);

	*vrot=out[1]/out[0];
	*sigR=sqrt(out[2]/out[0]); *sigp=sqrt(out[3]/out[0]); *sigz=sqrt(out[4]/out[0]); *sigRz=sqrt(out[5]/out[0]);
	return dens;
}


