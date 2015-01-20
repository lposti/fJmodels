#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "oct_int.h"
#include "potLEG.h"
#include "cuba.h"

#define LMAX 3
#define MAX(A,B) ((A)>(B)?(A):(B))

double df(double*,double*);
static double vscale,vmax;

double fL(double L){
	return 0;//tanh(L);
}
double rhoi(double *Rz,double VR,double Vphi,double Vz){//integrand for density etc
	double v[3]={VR,Vphi,Vz};
	return df(Rz,v);
}

double sinh_vgrid(double v){
	return vscale*sinh(asinh(vmax/vscale)*v/vmax);
}
/* to compute jacobian */
double cosh_vgrid(double v){
	return vscale*cosh(asinh(vmax/vscale)*v/vmax);
}
/*
 *  change of variables: (VR,Vphi,Vz) |---> (sinh(VR),sinh(Vphi),sinh(Vz))
 */
double rhoi_sinh(double *Rz,double VR,double Vphi,double Vz){
	double v[3]={sinh_vgrid(VR),sinh_vgrid(Vphi),sinh_vgrid(Vz)};
	double jacob=fabs(pow(asinh(vmax/vscale)/vmax,3)*cosh_vgrid(VR)*cosh_vgrid(Vphi)*cosh_vgrid(Vz));
	return jacob*df(Rz,v);
}
double rho_sinh(double R,double z,double *Vrot,double *sigmaR,double *sigmaP,double *sigmaZ,double *sigmaRZ,double sigma){
	double Phigl=Phi(R,z);
	double Rz[4]={R,z,Phigl,Phi(0,0)};
	// maximum speed = escape speed; scale speed = sigma
	double Ve=sqrt(-2*(Phigl-Phi(100,100))), Vze=Ve;
	vscale=sigma; vmax=Ve;
	double I0=oct_int(Rz,&rhoi_sinh,&fL,sinh_vgrid(-Ve),sinh_vgrid(Ve),sinh_vgrid(1e-6*Ve),sinh_vgrid(Ve),sinh_vgrid(0.),sinh_vgrid(Vze),
					  LMAX,Vrot,sigmaR,sigmaP,sigmaZ,sigmaRZ);
	*Vrot/=I0;// *sigmaP=sqrt(MAX(0,*sigmaP/I0-pow(*Vrot,2)));
	*sigmaP=sqrt(*sigmaP/I0);
	*sigmaR=sqrt(*sigmaR/I0); *sigmaZ=sqrt(*sigmaZ/I0); *sigmaRZ/=I0;
	return 4*I0;
}

/*
 * 19/01/15: Cuba integration
 */
struct CubaPars {double * Rz; double Ve;};
static int CubaInteg(const int *ndim, const double x[],
		  const int *ncomp, double f[], void * userdata){
	struct CubaPars * CP = (struct CubaPars *) (userdata);
	double * Rz = (CP->Rz);
	double Ve   = (CP->Ve);

	/* rescaling the integral to [0,1]x[0,1]x[0,1] */
	double V[3] = {1e-6+Ve*x[0], 1e-6-Ve+2.*Ve*x[1], Ve*x[2]};
	double jacob = 2.*Ve*Ve*Ve;

	//printf("%f %f %f %f %f\n",df(Rz,V),V[0],V[1],V[2],Ve);
	f[0] = df(Rz,V)*jacob;
	f[1] = f[0]*V[0]*V[0];
	f[2] = f[0]*V[1]*V[1];
	f[3] = f[0]*V[2]*V[2];

	return 0;
}
double rhoCuba(double R, double z, double * sigmaR, double * sigmaP, double * sigmaZ){
	double Phigl=Phi(R,z);
	double Rz[4]={R,z,Phigl,Phi(0,0)},Ve=sqrt(-2*(Phigl-Phi(100,100)));
	//return 4*oct_int(Rz,&rhoi,-Ve,Ve,1e-6*Ve,Ve,0,Vze,LMAX);

	struct CubaPars CP = {Rz,Ve};
	void * userdata = (void *) &CP;

	/* inputs */
	int NDIM=3,NCOMP=4,NVEC=4,
		VERBOSE=0,SEED=0,MINEVAL=0,MAXEVAL=50000,
		NSTART=1000,NINCREASE=500,NBATCH=1000,
		GRIDNO=0;
	double epsabs=1e-10,epsrel=1e-8;

	/* outputs */
	int neval,fail;
	double integral[NCOMP],error[NCOMP],prob[NCOMP];


	 Vegas(NDIM, NCOMP, CubaInteg, userdata, NVEC,
	    epsrel, epsabs, VERBOSE, SEED,
	    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
	    GRIDNO, NULL, NULL,
	    &neval, &fail, integral, error, prob);

	 /*
	  * Suave is alternative, but has been tweaked.
	  * Cuhre is not working properly

	int nreg;
	Suave(NDIM, NCOMP, CubaInteg, userdata, NVEC,
		 	    epsrel, epsabs, VERBOSE, SEED,
		 	    MINEVAL, MAXEVAL, 50000, 500, 4.,
		 	    NULL, NULL,
		 	    &nreg, &neval, &fail, integral, error, prob);

	 Cuhre(NDIM, NCOMP, CubaInteg, userdata, NVEC,
	 	    epsrel, epsabs, VERBOSE,
	 	    MINEVAL, MAXEVAL, 9,
	 	    NULL, NULL,
	 	    &nreg, &neval, &fail, integral, error, prob);
	*/

	 //double rho(double,double);
	 //printf("%e %e\n",4*integral[0],rho(R,z)); exit(1);
	 *sigmaR=sqrt(integral[1]/integral[0]);
	 *sigmaP=sqrt(integral[2]/integral[0]);
	 *sigmaZ=sqrt(integral[3]/integral[0]);
	 return 4*integral[0];
}

/*
 *  rho(R,z,*Vrot,*sigmaR,*sigmaP,*sigmaZ,*sigmaRZ,Vlim):
 *  adds the control over the integration limits in velocity space
 */
double rho(double R,double z,double *Vrot,double *sigmaR,double *sigmaP,double *sigmaZ,double *sigmaRZ,double Vlim){
	double Phigl=Phi(R,z);
	double Rz[4]={R,z,Phigl,Phi(0,0)},Ve=sqrt(-2*(Phigl-Phi(100,100))), Vze=Ve;
	double I0=oct_int(Rz,&rhoi,&fL,-Ve,Ve,1e-6*Ve,Ve,0,Vze,LMAX,Vrot,sigmaR,sigmaP,sigmaZ,sigmaRZ);
	*Vrot/=I0;// *sigmaP=sqrt(MAX(0,*sigmaP/I0-pow(*Vrot,2)));
	*sigmaP=sqrt(*sigmaP/I0);
	*sigmaR=sqrt(*sigmaR/I0); *sigmaZ=sqrt(*sigmaZ/I0); *sigmaRZ/=I0;
	return 4*I0;
}
double rho(double R,double z,double *Vrot,double *sigmaR,double *sigmaP,double *sigmaZ,double *sigmaRZ){
	// returns density, Vrot and all 3 dispersions
	double Phigl=Phi(R,z);
	double Rz[4]={R,z,Phigl,Phi(0,0)},Ve=sqrt(-2*(Phigl-Phi(100,100))), Vze=Ve;
	double I0=oct_int(Rz,&rhoi,&fL,-Ve,Ve,1e-6*Ve,Ve,0,Vze,LMAX,Vrot,sigmaR,sigmaP,sigmaZ,sigmaRZ);
	*Vrot/=I0;// *sigmaP=sqrt(MAX(0,*sigmaP/I0-pow(*Vrot,2)));
	*sigmaP=sqrt(*sigmaP/I0);
	*sigmaR=sqrt(*sigmaR/I0); *sigmaZ=sqrt(*sigmaZ/I0); *sigmaRZ/=I0;
	return 4*I0;
}
double rho(double R,double z,double *Vrot,double *sigmaR,double *sigmaP,double *sigmaZ){
	// returns density, Vrot and all 3 dispersions
	double Phigl=Phi(R,z);
	double Rz[4]={R,z,Phigl,Phi(0,0)},Ve=sqrt(-2*(Phigl-Phi(100,100))), Vze=Ve;
	double I0=oct_int(Rz,&rhoi,-Ve,Ve,1e-6*Ve,Ve,0,Vze,LMAX,Vrot,sigmaR,sigmaP,sigmaZ);
	*Vrot/=I0;// *sigmaP=sqrt(MAX(0,*sigmaP/I0-pow(*Vrot,2)));
	*sigmaP=sqrt(*sigmaP/I0);
	*sigmaR=sqrt(*sigmaR/I0); *sigmaZ=sqrt(*sigmaZ/I0);
	return 4*I0;
}
double rho(double R,double z,double *Vrot){// returns density and Vrot
	double Phigl=Phi(R,z);
	double Rz[4]={R,z,Phigl,Phi(0,0)},Ve=sqrt(-2*(Phigl-Phi(100,100))), Vze=Ve;
	double I0=oct_int(Rz,&rhoi,-Ve,Ve,1e-6*Ve,Ve,0,Vze,LMAX,Vrot);
	*Vrot/=I0;
	return 4*I0;
}
double rho(double R,double z){// returns density
	double Phigl=Phi(R,z);
	double Rz[4]={R,z,Phigl,Phi(0,0)},Ve=sqrt(-2*(Phigl-Phi(100,100))), Vze=Ve;
	return 4*oct_int(Rz,&rhoi,-Ve,Ve,1e-6*Ve,Ve,0,Vze,LMAX);
}
double nUi(double *Rz,double Vphi,double Vz){//integrand for U distrib
	double v[3]={Rz[4],Vphi,Vz};
	return df(Rz,v);
}
double nU(double R,double z,double VR){// returns U distribution
	double Phigl=Phi(R,z);
	double Rz[5]={R,z,Phigl,Phi(0,0),VR},Ve=sqrt(-2*(Phigl-Phi(100,100)));
	return oct_int(Rz,&nUi,.01*Ve,.7*Ve,-.7*Ve,.7*Ve,LMAX);
}
double nVi(double *Rz,double VR,double Vz){//integrand for V distrib
	double v[3]={VR,Rz[4],Vz};	
	return df(Rz,v);
}
double nV(double R,double z,double Vphi){// returns V distribution
	double Phigl=Phi(R,z);
	double Rz[5]={R,z,Phigl,Phi(0,0),Vphi},Ve=sqrt(-2*(Phigl-Phi(100,100)));
	return 2*oct_int(Rz,&nVi,-.7*Ve,.7*Ve,0,.7*Ve,LMAX);
}
double nWi(double *Rz,double VR,double Vphi){//integrand for W distrib
	double v[3]={VR,Vphi,Rz[4]};	
	return df(Rz,v);
}
double nW(double R,double z,double Vz){// returns W distribution
	double Phigl=Phi(R,z);
	double Rz[5]={R,z,Phigl,Phi(0,0),Vz},Ve=sqrt(-2*(Phigl-Phi(100,100)));
	return 2*oct_int(Rz,&nWi,-.7*Ve,.7*Ve,.01*Ve,.7*Ve,LMAX);
}
double nW(double R,double z,double Vz,double *Vrot){//returns W distribution and mean Vphi
	double Phigl=Phi(R,z);
	double Rz[5]={R,z,Phigl,Phi(0,0),Vz},Ve=sqrt(-2*(Phigl-Phi(100,100)));
	double I0=oct_int(Rz,&nWi,-.7*Ve,.7*Ve,.01*Ve,.7*Ve,LMAX,Vrot);
	*Vrot /=I0;
	return 2*I0;
}
