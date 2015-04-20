/*
 * Potential.cpp
 *
 *  Created on: Feb 15, 2015
 *      Author: L. Posti
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include "Potential.h"
#include "models.h"
#include "Grid.h"
#include "Utils.h"
#include "UtilsLeg.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"

Potential::Potential(const unsigned comp_in) {

	comp=comp_in;

	philP = mat<double>(NR,NPOLY);  PrP = mat<double>(NR,NPOLY);    Pr2P = mat<double>(NR,NPOLY);
	rhlP = mat<double>(NR,NPOLY);

	poly = mat<double>(NPOLY,NGAUSS);
	I_int = mat<double>(NR,NPOLY);
	I_ext = mat<double>(NR,NPOLY);

	ci  = arr<double>(NGAUSS); wi = arr<double>(NGAUSS); si = arr<double>(NGAUSS);
	pol = arr<double>(NPOLY);
	this->rhoGuess = 0;

	if (phil[0][0]<0) canEv = true;
	else canEv = false;

	// initializes Legendre Polynomials
	this->initLeg();
}

/*
 *  Selects the guess density from a list of classical
 *  methods.
 *
 *  INPUT:  string, can be ="Isochrone","Hernquist","NFW","NFWext"
 *  OUTPUT: void, initializes the private data rho.
 */
void Potential::selectGuessRho(const std::string rhoName){

	if      (rhoName=="Isochrone") this->rhoGuess = &rhoIsoch;
	else if (rhoName=="Hernquist") this->rhoGuess = &rhoHern;
	else if (rhoName=="NFW")       this->rhoGuess = &rhoNFW;
	else if (rhoName=="NFWext")    this->rhoGuess = &rhoNFWext;
	else {
		printf("\n\tCall selectGuessRho method with 'Isochrone' \n");
		exit(1);
	}
}

// rhl[n][i] now contains l=2*i Legendre coefficient of angular
// distribution of density on shell radius ar[n]
void Potential::computeGuessRhl(){

	for (int nr=0; nr<NR; nr++)
		for (int np=0; np<NPOLY; np++)
			rhlP[nr][np]=0.;

	for (int nr=0; nr<NR; nr++)
		for (int np=0; np<NPOLY; np++)
			for (int ng=0; ng<NGAUSS; ng++)
				rhlP[nr][np] += this->poly[np][ng] * this->rhoGuess(ar[nr]*this->si[ng],ar[nr]*this->ci[ng]);

	rhlP[0][0] = 2*this->rhoGuess(ar[0]*this->si[0],ar[0]*this->ci[0]);

	if (comp==1)
		for (int nr=0; nr<NR; nr++)
			for (int np=0; np<NPOLY; np++)
				rhl[nr][np]=rhlP[nr][np];
}

struct I_par { int l; double *rhl_h;};
double I_internal(double x, void * params){
	struct I_par * p = (struct I_par *) params;
	int l = p->l;
	double * rhl_h = p->rhl_h;

	return pow(x,l+2)*linterp<double> (ar,rhl_h,NR,x) ;
}
double I_external(double x, void * params){
	struct I_par * p = (struct I_par *) params;
	int l = p->l;
	double * rhl_h = p->rhl_h;

	// after the last grid point we assume Keplerian regime: rho~r^-4
	double rhol;
	if (x>ar[NR-1]) rhol = MAX (0., rhl_h[NR-1]*pow(ar[NR-1]/x,4));
	else rhol = linterp<double> (ar,rhl_h,NR,x);

	return rhol/pow(x,l-1);
}

/*
 *  Compute internal and external integrals to get
 *  the coefficients of the Potential multipole expansion
 */
#define EPSABS   1e-4
#define EPSREL   1e-3
#define WORKSIZE 1000
void Potential::computeInts(double **rhlH){

	gsl_integration_workspace * w
		= gsl_integration_workspace_alloc (WORKSIZE);

	gsl_function F,G;
	F.function = &I_internal;
	G.function = &I_external;

	double res,err;
	for (int np=0; np<NPOLY; np++){

		// array to be used by the linear interpolation scheme
		double *rhl_h = arr<double> (NR);

		for (int i=0; i<NR; i++) rhl_h[i] = rhlH[i][np];

		for (int nr=0; nr<NR; nr++){

			struct I_par par = {2*np, rhl_h};
			F.params = &par; G.params = &par;

			// Internal integral
			gsl_integration_qags (&F, 0., ar[nr], EPSABS, EPSREL, WORKSIZE, w, &res, &err);
			this->I_int[nr][np] = res;

			// External Integral
			gsl_set_error_handler_off();
			double epsabs=EPSABS,epsrel=EPSREL;
			int status=1;

			while(status!=GSL_SUCCESS){
				status = gsl_integration_qagiu (&G, ar[nr], epsabs, epsrel, WORKSIZE, w, &res, &err);
				epsabs*=2; epsrel*=2;
			}
			this->I_ext[nr][np] = res;

		}
		delArr(rhl_h);
	}

	gsl_integration_workspace_free (w);
}

/*
 *  Compute the coefficients of the asymptotic expansion of the potential
 *
 *  INPUT:  none, calls the function computeInts() if not already done
 *  OUTPUT: none, fills the global arrays phil, Pr, Pr2
 */
void Potential::computePhil(double **rhlH){

	// Compute Internal and External Integrals, if not done before.
	this->computeInts(rhlH);

	for (int np=0; np<NPOLY; np++){
		int l=2*np; double G=1;
		for(int n=0; n<NR; n++){
			philP[n][np]=-TPI*G*(this->I_int[n][np]/pow(ar[n],l+1)+
						(this->I_ext[n][np])*pow(ar[n],l));
			PrP[n][np]=-TPI*G*(-(l+1)*this->I_int[n][np]/pow(ar[n],l+2)+
					  l*(this->I_ext[n][np])*pow(ar[n],l-1));

			Pr2P[n][np]=-TPI*G*((l+2)*(l+1)*this->I_int[n][np]/pow(ar[n],l+3)+
				   (l-1)*l*(this->I_ext[n][np])*pow(ar[n],l-2)
				   -(2*l+1)*rhlH[n][np]);
			}
	}

	if (comp==1)
		for (int n=0; n<NR; n++)
			for (int np=0; np<NPOLY; np++){
				phil[n][np]=philP[n][np];
				Pr[n][np]  =PrP[n][np];
				Pr2[n][np] =Pr2P[n][np];
			}
	/*
	else
		for(int n=0; n<NR; n++)
			for (int np=0; np<NPOLY; np++){
				phil[n][np]+=philP[n][np];
				Pr[n][np]  +=PrP[n][np];
				Pr2[n][np] +=Pr2P[n][np];
			}
	*/

	canEv = true; // automatically reset the potential to re-compute the integrals
}

/*
 *  Evaluate potential at cylindrical (R,z)
 *  Needs Legendre coefficients phil to be present
 */
double Potential::operator() (double R,double z){

	// Compute phils, if never done before.
	if (!this->canEv) this->computePhil();

	double pol[NPOLY],phip[NPOLY];
	double c=z/sqrt(R*R+z*z), r=sqrt(R*R+z*z);

	intpo2<double>(r,phip); evenlegend<double>(pol,c);

	double phi=phip[0];
	for(int np=1; np<NPOLY; np++) phi+=phip[np]*pol[np];
	return phi;
}

/*
 *  Spherical derivatives of the Potential
 */
double Potential::dr(const double r,const double theta,double *d2Phidr2){
	double pol[NPOLY],phip[NPOLY],dphip[NPOLY],d2phip[NPOLY];
	double grad1=0,grad2=0;
	evenlegend(pol,cos(theta));
	intpo2(r,phip,dphip,d2phip);
	for(int np=0; np<NPOLY; np++){
		grad1+=dphip[np]*pol[np];
		grad2+=d2phip[np]*pol[np];
	}
	(*d2Phidr2)=grad2; return grad1;
}
double Potential::dtheta(const double r, const double c) {
    //returns the first derivative of the whole potential wrt spherical theta
	double dphitheta = 0;
	double dpol[NPOLY], phinterp[NPOLY];
	intpo2(r, phinterp);
	for(int k=0; k<NPOLY; k++) dpol[k] = dlegend<double>(c, 2*k); //fills the array dpol with derivatives of the even legendre polys
	for(int j=0; j<NPOLY; j++) dphitheta += phinterp[j]*dpol[j]; //this is the derivative wrt cos(theta), need to then multiply by -sin(theta)
	dphitheta *= -sqrt(1 - c*c); //im taking the positive root here - any issues with this?
	return dphitheta;
}

/*
 *  Cylindrical derivatives of the Potential
 */
double Potential::dR(const double R, const double z) {//derivative of Phi with respect to cylindrical R
	double r = sqrt(z*z + R*R);
	if (r<ar[0]) r=ar[0]; // 25\02\15 LP edit: avoiding nans...
	double s=R/r, c=z/r;
	if(c==1) c-=1e-8; // 15\05\14 LP edit: just removing tiny value for z>>r (otherwise nans)
	double theta=acos(c), d2;
	return this->dtheta(r, c)*(c/r) + this->dr(r, theta, &d2)*s; //chain rule
}
double Potential::dz(const double R, const double z) {//derivative of Phi with respect to cylindrical z
	double r = sqrt(z*z + R*R);
	if (r<ar[0]) r=ar[0]; // 25\02\15 LP edit: avoiding nans...
	double s=R/r, c=z/r;
	if(c==1) c-=1e-8; // 15\05\14 LP edit: just removing tiny value for z<<r (otherwise nans)
	double theta=acos(c), d2;
	return this->dtheta(r, c)*(-s/r) + this->dr(r, theta, &d2)*c; //chain rule
}


/*
 *  Fills the arrays needed for the multipole expansion
 */
void Potential::initLeg(){
	gauleg<double>(0,1,ci,wi);
	fillPoly<double>(this->si,this->ci,this->wi,this->pol,this->poly);
}
