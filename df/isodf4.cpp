#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uv_orb3.h"
#include "tables3.h"
#include "press.h"
#include "ELz_grid.hh"
//#include "eddington.hh"
#include "moments.h"
#include "leg_pot3.h"
#include "ini_potLEG.h"
#include "oct_int.h"
#include "isodf4.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_interp.h"

#define TPI 6.283185307179586
#define SQRT2 1.4142135623730951
#define SMALLEST 10		// see tables3.cpp
#define DFGRID 200

static const double GM=1,a=1;
static double dr,dphi,dz;
#ifdef HJGJ
static double dr_g,dphi_g,dz_g;
#endif
static int count;

void setdf(double dr_in, double dphi_in,double dz_in){
	dr=dr_in; dphi=dphi_in; dz=dz_in;
}

#ifdef HJGJ
void setdf(double dr_in, double dphi_in, double dz_in,
		double dr_g_in, double dphi_g_in, double dz_g_in){
	dr=dr_in;     dphi=dphi_in;     dz=dz_in;
	dr_g=dr_g_in; dphi_g=dphi_g_in; dz_g=dz_g_in;
}
#endif

int myisnan(double x){if (x==x) return 0; else return 1;}

double hj_COmRatio(double Jr, double Jphi, double Jz){
	return dr*Jr+dphi*Jphi+dz*Jz;
}

#ifdef HJGJ
double gj_COmRatio(double Jr, double Jphi, double Jz){
	return dr_g*Jr+dphi_g*Jphi+dz_g*Jz;
}
#endif

double hj(double Jr, double Lz, double Jz, double R, double J0){
		double kappa,nu,Omegac;
		double Rc,Lc,Lc_old=0.;
		count++;

		ELz_grid grid;
		Lc=Jr+sqrt(2.)*Lz+sqrt(2.)*Jz;

		//printf("\n>> R=%f, %e %e %e Omegas %e %e\n",R,Jr,Lz,Jz,uvorb.Omegau,uvorb.Omegav);
		int nn=0;
		while(fabs(Lc-Lc_old)>1e-3){
			Lc_old=Lc;

			Rc=grid.GetRc_ini(Lc*Lc,R);
			getfreqs_ini(Rc,&kappa,&nu,&Omegac);
			Lc=Jr+Omegac/kappa*Lz+nu/kappa*Jz;
			if (myisnan(Lc)==1) printf("Lc is NANANANANANAN\n");
			nn++; if (nn>50) printf("Difficulty in computing h(J)... Lc=%f Lc_old=%f\n",Lc,Lc_old);
		}
		//fprintf(fdf,"%lf %lf %lf %lf %lf\n",Lz,pow(Lc,-2),R,nu/kappa,Omegac/kappa);


		if (kappa==0.) return -1;
		//return (Jr+J0)+.4*Omegac/kappa*(Lz+J0)+1.2*nu/kappa*(Jz+J0);;
		//	spherical
		//return Lc=Jr+J0+Omegac/kappa*(Lz+J0)+nu/kappa*(Jz+J0);
		//return Jr+fabs(Lz)+Jz;

		//	azimuthal bias deltaphi=.4,deltaz=1.2
		//  radial bias deltaphi=1.4,deltaz=1.5
		double deltar,deltaphi=.4,deltaz=1.2;
		//deltar=1.-Omegac/kappa*(deltaphi-Omegac/kappa)-nu/kappa*(deltaz-nu/kappa);
		deltar=1.-pow(Omegac/kappa,1)*(deltaphi-1.)-pow(nu/kappa,1)*(deltaz-1.);

		Lc=deltar*Jr+deltaphi*Lz+deltaz*Jz;
		//return Lc;

		//  modification for low ellipticities
		//double deltazero=(pow(Omegac/kappa,2)+pow(nu/kappa,2)-deltaz*nu/kappa+1.)/(Omegac/kappa+1.);
		//double deltazero=1.-pow(nu,2)/(pow(Omegac,2)+pow(kappa,2))*(deltaz-1.);
		double deltazero=1-(deltaz-1)/(1+kappa/Omegac);
		double psi=tanh(Lc/sqrt(GM*a));
		deltaphi=(1.-psi)*deltazero+psi*deltaphi;
		deltar=(1.-psi)*deltazero+psi*deltar;

		// use for spherical
		//deltaz=deltaphi;

		Lc=deltar*(Jr+J0)+deltaphi*Omegac/kappa*(Lz+J0)+deltaz*nu/kappa*(Jz+J0);
		//*Om=Omegac/kappa;

		return Lc;
}
/*
 *  f(J)=f[H(J)] \propto h(J)^-2, where h(J)=Jr+Omega_phi/Omega_r*Jphi+Omega_z/Omega_r*z
 *  Isothermal model
 */
double hj_isoth(double Jr, double Lz, double Jz, double R){
	double mass=1e2, J0=1, JM=21.15;

#ifdef HJGJ
	double hJ=hj_COmRatio(Jr,Lz,Jz),gJ=gj_COmRatio(Jr,Lz,Jz);
	return MAX(0.,1./(pow(J0+gJ,2.))/mass  -
					1./(pow(J0+3.*JM,2.))/mass);
#else
#	ifdef CONSTOMRATIO
		double J=hj_COmRatio(Jr,Lz,Jz);
#	else
		double J=hj(Jr,Lz,Jz,R,J0);
#	endif /* not CONSTOMRATIO */
		return MAX(0.,pow(J+J0,-2.5)/mass-pow(J0+JM,-2.5)/mass);
#endif /* not HJGJ */
}

/*
 *  f(J)=f[H(J)] \propto h(J)^-5
 *  Isochrone model
 */
double hj_isoch(double Jr, double Lz, double Jz, double R){
	double mass=1.4e1,J0=1,JM=11.15;

#ifdef HJGJ
	double hJ=hj_COmRatio(Jr,Lz,Jz),gJ=gj_COmRatio(Jr,Lz,Jz);
	return MAX(0.,1./(pow(J0+gJ,5.))/mass  -
				1./(pow(J0+3.*JM,5.))/mass);
#else
#	ifdef CONSTOMRATIO
		double J=hj_COmRatio(Jr,Lz,Jz); //Jr+.5*(Lz+4.*Jz);
#	else
		double J=hj(Jr,Lz,Jz,R,J0);
#	endif /* not CONSTOMRATIO */
		return MAX(0.,pow(J+J0,-5)/mass-pow(J0+3*JM,-5)/mass);
#endif /* not HJGJ */
}

/*
 *  f(J)=f[H(J)] \propto h(J)^-5/3 * h(J)^-5+5/3
 *  Hernquist model
 */
double hj_hernq(double Jr, double Lz, double Jz, double R){
	double mass=0.3e2,J0=1,JM=11.15;

#ifdef HJGJ
	double hJ=hj_COmRatio(Jr,Lz,Jz),gJ=gj_COmRatio(Jr,Lz,Jz);
	return MAX(0.,pow(J0+hJ,5./3.)/(pow(hJ,5./3.)*pow(J0+gJ,5.-5./3.))/mass  -
				pow(J0+3.*JM,5./3.)/(pow(3.*JM,5./3.)*pow(J0+3.*JM,5.-5./3.))/mass);
#else
#	ifdef CONSTOMRATIO
		double J=hj_COmRatio(Jr,Lz,Jz);
#	else
		double J=hj(Jr,Lz,Jz,R,J0);
#	endif /* not CONSTOMRATIO */
	return MAX(0.,pow(J,-5./3.)*pow(J0+J,-5+5./3.)/mass-pow(JM,-5./3.)*pow(J0+JM,-5+5./3.)/mass);
#endif /* not HJGJ */

}

/*
 *  f(J)=f[H(J)] \propto h(J)^-5/3 * h(J)^-5+5/3
 *  NFW model
 */
double hj_nfw(double Jr, double Lz, double Jz, double R){
	double mass=1.4e1,J0=1,JM=11.15;

#ifdef HJGJ
	double hJ=hj_COmRatio(Jr,Lz,Jz),gJ=gj_COmRatio(Jr,Lz,Jz);
	return MAX(0.,pow(J0+hJ,5./3.)/(pow(hJ,5./3.)*pow(J0+gJ,3.-5./3.))/mass  -
				pow(J0+3.*JM,5./3.)/(pow(3.*JM,5./3.)*pow(J0+3.*JM,3.-5./3.))/mass);
#else
#	ifdef CONSTOMRATIO
		double J=hj_COmRatio(Jr,Lz,Jz);
#	else
		double J=hj(Jr,Lz,Jz,R,J0);
#	endif /* not CONSTOMRATIO */
	return MAX(0.,pow(J,-5./3.)*pow(J0+J,-3+5./3.)/mass-pow(JM,-5./3.)*pow(J0+JM,-3+5./3.)/mass);
#endif /* not HJGJ */
}

double df(double *x,double *v){
	double R=x[0],Lz=R*v[1],Phigl=x[2];

	/* Hamiltonian: here is needed to compute Delta */
	double H=.5*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2))+Phigl;
	if (H>=0) return 0.;

	/* get the actions */
	double p[2]={v[0],v[2]},Jr,Jz;
	uv_orb uvorb(Deltafn(H),Lz,Phigl,x,p);
	/*
	 * INTERPOLATION SCHEME NOT USED ANYMORE!
	 *
	 * int nt=int_acts(Lz,uvorb.E,uvorb.Er,&Jr,&Jz);
	 * if(nt==0){
	 * 	Jr=uvorb.Ju(); Jz=uvorb.Jv();// Jz=JzI;
	 * }
	 */

	Jr=uvorb.Ju(); Jz=uvorb.Jv();

	//uvorb.GetFreqs();
	//double Omegar=uvorb.Omegau,Omegaphi=uvorb.Omegaphi;

#if defined(HERNQUIST) || defined(NFW)			// use the same procedure for Hernquist and NFW models
	return hj_hernq(Jr,Lz,Jz,R);

#elif defined ISOTHERMAL
	return hj_isoth(Jr,Lz,Jz,R);

#elif defined ISOCHRONE
	return hj_isoch(Jr,Lz,Jz,R);

#endif
}
