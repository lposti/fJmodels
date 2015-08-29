/*
 * df.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: L. Posti
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

void setDF(struct fJParams * fJP, Potential *p, const unsigned comp){

	if (comp==1){

		if (fJP->modName=="Hernquist")				{alpha[0] = fJP->alpha = 1.667;   alpha[1] = fJP->beta = 5.;}
		else if (fJP->modName=="Isochrone") 		{alpha[0] = fJP->alpha = 0.0;     alpha[1] = fJP->beta = 5.;}
		else if (fJP->modName=="NFW")				{alpha[0] = fJP->alpha = 1.667;   alpha[1] = fJP->beta = 3.;}
		else if (fJP->alpha != UNSET and fJP->beta != UNSET)     {alpha[0] = fJP->alpha;   alpha[1] = fJP->beta;}
		else {printf("\n--> ERROR: either set parameter 'model' or 'alpha' and 'beta'!\n"); exit(1);}

		J0=sqrt(fJP->mass*fJP->r0); mass=fJP->mass; chi=fJP->chi;
		setDF(fJP->dphi_h_in,fJP->dz_h_in,fJP->dphi_g_in,fJP->dz_g_in,p);

	} else if (comp==2){

		if (fJP->modName2=="Hernquist")				{alpha[0] = fJP->alpha_2 = 1.667;   alpha[1] = fJP->beta_2 = 5.;}
		else if (fJP->modName2=="Isochrone")		{alpha[0] = fJP->alpha_2 = 0.0;     alpha[1] = fJP->beta_2 = 5.;}
		else if (fJP->modName2=="NFW")				{alpha[0] = fJP->alpha_2 = 1.667;   alpha[1] = fJP->beta_2 = 3.;}
		else if (fJP->alpha != UNSET and fJP->beta != UNSET)     {alpha[0] = fJP->alpha_2;   alpha[1] = fJP->beta_2;}
		else {printf("\n--> ERROR: either set parameter '2:model' or '2:alpha' and '2:beta'!\n"); exit(1);}

		J0=sqrt(fJP->mass_2*fJP->r0_2); mass=fJP->mass_2; chi=fJP->chi_2;
		setDF(fJP->dphi_h_in2,fJP->dz_h_in2,fJP->dphi_g_in2,fJP->dz_g_in2,p);

	}

	printf("===========  DF Set as: ==============\n a=%f b=%f \n",alpha[0],alpha[1]);
}

inline void hJgJ(double * hg, const double Jr, const double Jphi, const double Jz){
	hg[0]=Jr+dphi_h*fabs(Jphi)+dz_h*Jz; hg[1]=Jr+dphi_g*fabs(Jphi)+dz_g*Jz;
}

/*
 *  f(J)=f[h(J),g(J)]
 */
double df_hg(const double Jr, const double Jphi, const double Jz){

	double hg[2]; hJgJ(hg,Jr,Jphi,Jz);
	return mass/pow(J0,3) * MAX(0., pow(1.+J0/hg[0],alpha[0])/(pow(1.+hg[1]/J0,alpha[1]))) / (pow(TPI,3.));
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

	double DFeven = df_hg(Jr,Jphi,Jz);
	double k=K_SATOH;
	//double rc=RcE<double>(H, 1., gP);
	//if ((1.-k)*DFeven + k * tanh(chi*Jphi/J0)*DFeven>1e-6)
	//	fprintf(f_lz, "%e %e %e %e %e %e \n", R, log10((1.-k)*DFeven + k * tanh(chi*Jphi/J0)*DFeven), H, Jphi, rc*sqrt(rc*gP->dR(rc, 0.)), Lz / (rc*sqrt(rc*gP->dR(rc, 0.))));

	if (chi==0.) return DFeven;
	else if (chi>0.){
		//double k=0.5;  // 0.625 for ngc6427
		return (1.-k)*DFeven + k * tanh(chi * Jphi / J0)*DFeven;
	}
	else {
		printf("\n  --ERROR: chi must be >= 0!");
		exit(1);
	}
}
