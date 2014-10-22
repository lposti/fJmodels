#include <math.h>
#include <stdio.h>
#include "oct_int.h"
#include "potLEG.h"

#define LMAX 3
#define MAX(A,B) ((A)>(B)?(A):(B))

double df(double*,double*);
static double vscale,vmax;

double fL(double L){
	return tanh(L);
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