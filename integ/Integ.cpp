/*
 * integ.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: L. Posti
 */

#include <math.h>
#include "Grid.h"
#include "DF.h"
#include "Potential.h"
#include "GaussQuad.h"
#include "Utils.h"
#include "UtilsLeg.h"


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

/* Edge-on density projection */
struct SigmaIntegPar {double R; double z;};
double SigmaInteg (const double y, void * params){

	struct SigmaIntegPar * par = (struct SigmaIntegPar*) params;
	double Rm = par->R;
	double z = par->z;

	double r = Rm + y * (ar[NR-1]-Rm);
	double jacob = (ar[NR-1]-Rm);

	return jacob*ev_dens(r, z);
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

struct LPPar { double * Rz; double Ve; Potential *p; double vx;};
double lineProfile(const double vy, const double x, const double vz, void * params){

	struct LPPar * par = (struct LPPar*) params;
	double Ve = par->Ve;
	double *Rz = par->Rz;
	Potential *p = par->p;
	double vx = par->vx;

	double R = sqrt(pow(Rz[0],2) + pow(x * (ar[NR-1]),2)), z = Rz[1], y = Rz[0];
	// double vphi = sgn<double>(vx)*sqrt(vx*vx+vy*vy) * sin(atan(y/x));
	// double vr = vy;
	double v_abs = sqrt(vx*vx + vy*vy);
	double alpha = acos(vx/v_abs), theta = atan2(y, x);
	double vr = v_abs * cos(alpha+theta), vphi = v_abs * sin(alpha+theta);

	//Rz[0]=R;
	//Rz[2]=(*p)(R,z);
	//Ve=sqrt(-2*((*p)(R,z)-(*p)(100 * ar[NR-1], 100 * ar[NR-1])));

	/* rescaling the integral to [0,1]x[-1,1]x[0,1] */
	double V[3] = {1e-6+Ve*vr, 1e-6+Ve*vphi, 1e-6+Ve*vz};
	double jacob = Ve*Ve*(fabs(x) * pow(ar[NR-1],2))/R;

	//if (x<0.)
	//	printf("%e %e %e %e %e %e\n",vr,vphi,vz,vx,vy,Ve);

	return jacob*df(Rz,V);
}

/*
 *  rho from DF integration
 */
double rhofDF(double R, double z, Potential *p){
	double Phi_h=(*p)(R,z);
	double Rz[3]={R,z,Phi_h},Ve=sqrt(-2*(Phi_h-(*p)(10 * ar[NR-1],10 * ar[NR-1])));

	struct IntegPar par = {&Rz[0],Ve};
	return 4*Int3D_011101(&rhoInteg,&par);
}

/*
 *  rho from DF integration: with sigma(R,z,phi,Rz) computed
 */
double rhofDF(double R, double z, Potential *p, double * vrot, double * sigR,
		      double * sigp, double * sigz, double * sigRz){
	double Phi_h=(*p)(R,z);
	double Rz[3]={R,z,Phi_h},Ve=sqrt(-2*(Phi_h-(*p)(10 * ar[NR-1],10 * ar[NR-1])));

	struct IntegPar par = {&Rz[0],Ve};
	double * out = arr<double>(6);
	//double dens=4*Int3D_011101_vec(&rhoInteg_vec,&par,out);
	double dens=Int3D_111111_vec(&rhoInteg_vec,&par,out);

	*vrot=out[1]/out[0];
	*sigR=sqrt(out[2]/out[0]); *sigp=sqrt(out[3]/out[0]); *sigz=sqrt(out[4]/out[0]); *sigRz=sqrt(fabs(out[5])/out[0]); // *sigRz=sqrt(MAX(0., out[5])/out[0]);
	return dens;
}

/*
 *  projected density: assuming edge-on line-of-sight
 */
double SigmaDF(double R, double z){

	struct SigmaIntegPar par = {R,z};
	return 2*Int1D_01(&SigmaInteg, &par);
}

/*
 *  line profile integration
 */
double line_profile(double R, double z, double vx, Potential *p){

	double Phi_h=(*p)(R,z);
	double Rz[3]={R,z,Phi_h},Ve=sqrt(-2*(Phi_h-(*p)(10 * ar[NR-1],10 * ar[NR-1])));

	double Sigma = SigmaDF(R,z);

	struct LPPar par = {&Rz[0], Ve, p, vx};
	return Int3D_111111(&lineProfile, &par) / Sigma;
}
