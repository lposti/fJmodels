/*
 * df.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: morpheus
 */

#include <cstring>
#include <math.h>
#include "DF.h"
#include "readParam.h"
#include "models.h"
#include "Utils.h"
#include "Tabulate.h"
#include "uvOrb.h"
#include "Potential.h"

static double alpha[2] = {5./3.,5.};
static double dphi_h, dz_h, dphi_g, dz_g;
static Potential *gP;
static std::string modClassic="Hernquist";

inline void setDF(const double dphi_h_in, const double dz_h_in,
                  const double dphi_g_in, const double dz_g_in){

	dphi_h = dphi_h_in; dz_h = dz_h_in;
	dphi_g = dphi_g_in; dz_g = dz_g_in;
}

void setDF(const double dphi_h_in, const double dz_h_in,
                  const double dphi_g_in, const double dz_g_in,
                  Potential *p){
	gP = p;
	setDF(dphi_h_in,dz_h_in,dphi_g_in,dz_g_in);
}

void setDF(struct fJParams fJP, Potential *p){

	modClassic = fJP.modName;
	if (modClassic=="Hernquist") 	  {alpha[0] = 5./3.; alpha[1]=5.;}
	else if (modClassic=="Isochrone") {alpha[0] = .0;    alpha[1]=5.;}
	else {printf("\n--> ERROR: currently only Isochrone and Hernquist models supported!\n"); exit(1);}

	setDF(fJP.dphi_h_in,fJP.dz_h_in,fJP.dphi_g_in,fJP.dz_g_in,p);
}

inline void hJgJ(double * hg, const double Jr, const double Jphi, const double Jz){
	hg[0]=Jr+dphi_h*fabs(Jphi)+dz_h*Jz; hg[1]=Jr+dphi_g*fabs(Jphi)+dz_g*Jz;
}

/*
 *  f(J)=f[h(J),g(J)]
 */
double df_hg(const double Jr, const double Jphi, const double Jz){

	double hg[2]; hJgJ(hg,Jr,Jphi,Jz);
	return mass/pow(J0,3)*MAX(0.,pow(1.+J0/hg[0],alpha[0])/(pow(1.+hg[1]/J0,alpha[1]))) / (pow(TPI,3.));
}

double df(const double *x, const double *v){

	double R=x[0],Lz=R*v[1],Phi_h=x[2];

	/* Hamiltonian: here is needed to compute Delta */
	double H=.5*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2))+Phi_h;
	if (H>=0) return 0.;

	/* get the actions */
	double p[2]={v[0],v[2]},Jr,Jphi,Jz;
	uvOrb uvorb(Deltafn(H),Lz,Phi_h,x,p,gP);
	Jr=uvorb.Ju(); Jz=uvorb.Jv(); Jphi=Lz;

	return df_hg(Jr,Jphi,Jz);
}
