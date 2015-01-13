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
#include "potLEG.h"
#include "isodf4.h"
#include "print.h"

#define TPI 6.283185307179586
#define SQRT2 1.4142135623730951
#define SMALLEST 10		// see tables3.cpp
#define DFGRID 200

static const double GM=1,a=1;
static double dr,dphi,dz;
double mass,J0,r0;
#ifdef HJGJ
static double dr_g,dphi_g,dz_g;
#endif
static int count;


void setdf(double dr_in, double dphi_in,double dz_in){
	dr=dr_in; dphi=dphi_in; dz=dz_in;
}

void setMJ0(double M_in, double J0_in){
	mass=M_in; J0=J0_in; r0=pow(J0_in,2)/M_in;
	// G=1 here
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
	return dr*Jr+dphi*fabs(Jphi)+dz*Jz;
}

#ifdef HJGJ
double gj_COmRatio(double Jr, double Jphi, double Jz){
	return dr_g*Jr+dphi_g*fabs(Jphi)+dz_g*Jz;
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

#ifdef HJGJ
	double gJ=gj_COmRatio(Jr,Lz,Jz);
	double rt=100,Jt=sqrt(rt*pow(TPI,3.));
	return mass/pow(J0,3)*MAX(0.,1./(pow(1+gJ/J0,2))-1./(pow(1.+Jt/J0,2))) / pow(TPI,3.);

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

#ifdef HJGJ
	double gJ=gj_COmRatio(Jr,Lz,Jz);
	return mass/pow(J0,3)*MAX(0.,1./(pow(J0+gJ,5.5))) / pow(TPI,3.);

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

#ifdef HJGJ
	double hJ=hj_COmRatio(Jr,Lz,Jz),gJ=gj_COmRatio(Jr,Lz,Jz);
	if (hJ<1e-4) return 0.;
	else return mass/pow(J0,3)*MAX(0.,pow(1.+J0/hJ,1.5)/(pow(1.+gJ/J0,5.))) / (pow(TPI,3.)) ;

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

#ifdef HJGJ
	double hJ=hj_COmRatio(Jr,Lz,Jz),gJ=gj_COmRatio(Jr,Lz,Jz);
	return mass/pow(J0,3)*MAX(0.,pow(1.+J0/hJ,5./3.)/(pow(1.+gJ/J0,2.9))) / pow(TPI,3.) ;

#else
#	ifdef CONSTOMRATIO
		double J=hj_COmRatio(Jr,Lz,Jz);
#	else
		double J=hj(Jr,Lz,Jz,R,J0);
#	endif /* not CONSTOMRATIO */
	return MAX(0.,pow(J,-5./3.)*pow(J0+J,-3+5./3.)/mass-pow(JM,-5./3.)*pow(J0+JM,-3+5./3.)/mass);
#endif /* not HJGJ */
}

/*
 *  f(J)=f[H(J)] \propto h(J)^-2 * h(J)^-5+2
 *  Jaffe model
 */
double hj_jaffe(double Jr, double Lz, double Jz, double R){

#ifdef HJGJ
	double hJ=hj_COmRatio(Jr,Lz,Jz),gJ=gj_COmRatio(Jr,Lz,Jz);
	return mass/pow(J0,3)*MAX(0.,pow(1.+J0/hJ,1.885)/(pow(1.+gJ/J0,5.))) / (pow(TPI,3.)) ;

#else
#	ifdef CONSTOMRATIO
		double J=hj_COmRatio(Jr,Lz,Jz);
#	else
		double J=hj(Jr,Lz,Jz,R,J0);
#	endif /* not CONSTOMRATIO */
	return MAX(0.,pow(J,-5./3.)*pow(J0+J,-5+5./3.)/mass-pow(JM,-5./3.)*pow(J0+JM,-5+5./3.)/mass);
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

#ifdef HERNQUIST			// use the same procedure for Hernquist and NFW models
	double DF = hj_hernq(Jr,Lz,Jz,R);

#elif defined ISOTHERMAL
	double DF = hj_isoth(Jr,Lz,Jz,R);

#elif defined ISOCHRONE
	double DF = hj_isoch(Jr,Lz,Jz,R);

#elif defined NFW
	double DF = hj_nfw(Jr,Lz,Jz,R);

#elif defined JAFFE
	double DF = hj_jaffe(Jr,Lz,Jz,R);
#endif

#ifdef PRINTDFH
	printDFH(Jr,Lz,Jz,DF,H);
#endif
#ifdef PRINTXVJ
	printXVJ(x,v,Jr,Lz,Jz);
#endif

	return DF;
}
