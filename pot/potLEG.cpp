#include <math.h>
#include "leg_pot2.h"

double Phi(double R,double z){
	return phileg(R,z);
}
double PhiR(double R){
	return phileg(R,0);
}
double dPhiR(double R, double z) {//derivative of Phi with respect to cylindrical R
	double r = sqrt(z*z + R*R), s=R/r, c=z/r;
	if(c==1) c-=1e-8; // 15\05\14 LP edit: just removing tiny value for z>>r (otherwise nans)
	double theta=acos(c), d2;
	return dPhitheta(r, c)*(c/r) + dPhir(r, theta, &d2)*s; //chain rule
}
double dPhiRnum(double R, double z) {return dPhiR(R,z);}
double dPhiz(double R, double z) {//derivative of Phi with respect to cylindrical z
	double r = sqrt(z*z + R*R), s=R/r, c=z/r;
	if(c==1) c-=1e-8; // 15\05\14 LP edit: just removing tiny value for z<<r (otherwise nans)
	double theta=acos(c), d2;
	return dPhitheta(r, c)*(-s/r) + dPhir(r, theta, &d2)*c; //chain rule
}
double dPhiR(double R){
	return dPhiR(R,0);
}
double dPhi(double *x2,double *f2){
	double R=x2[0],z=x2[1];
	f2[0]=dPhiR(R,z);
	f2[1]=dPhiz(R,z);
	return phileg(R,z);
}
double d2PhiR(double R) {
    //returns the second derivative of Phi wrt cylindrical R in the plane
	double r = R, theta=asin(1), d2Phidr2;
	dPhir(r, theta, &d2Phidr2);
	return d2Phidr2;
}
double d2Phiz(double R) {
    //returns the second derivative of Phi wrt cylindrical z
	double r = R, theta=asin(1), d2;
	return (1./r)*dPhir(r, theta, &d2) + (1./(r*r))*d2Phitheta(r, 0);
}
void getfreqs(double R,double *kappa,double *nu,double *Omega){
	(*kappa)=sqrt(d2PhiR(R)+3/R*dPhiR(R,0));
	(*nu)=sqrt(d2Phiz(R));
	(*Omega)=sqrt(dPhiR(R,0)/R);
}
