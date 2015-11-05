/*
 * models.cpp
 *
 *  Created on: Feb 15, 2015
 *      Author: L. Posti
 */

#include <math.h>
#include "models.h"
#include "Utils.h"
#include "readParam.h"

double mass,J0,r0,chi;
double q,q2;

void SetModel(struct fJParams fJP, const unsigned comp){

	if (comp==1){
		mass=fJP.mass; r0=fJP.r0;
		q=fJP.q; q2=q*q;
	} else if (comp==2) {
		mass=fJP.mass_2; r0=fJP.r0_2;
		q=fJP.q_2; q2=q*q;
	}
}

double logFlat(double R, double z){
	double v0=.15, Rc=.1;
	return v0*v0/(FPI*q2) * (Rc*Rc*(2*q2+1)+R*R+(2.-1./q2)*z*z) /
			pow(Rc*Rc+R*R+z*z/q2,2);
}

double MN(double R, double z){
	double b=.3, a=1, M=1;
	double b2=b*b, a2=a*a;
	return b2*M/FPI * (a*R*R + (a+3.*sqrt(z*z+b2))*pow(a+sqrt(z*z+b2),2) ) /
			pow(R*R + (a+sqrt(z*z+b2)), 5./2.) / pow(z*z+b2,3./2.);
}

double rhoHern(double R, double z){
	//return logFlat(R,z);
	//return MN(R,z);
	double m=sqrt(R*R+z*z/q2);
	return mass/(m/r0*pow(r0+m,3))/TPI;
}

double rhoIsoch(double R, double z){
	double m=sqrt(R*R+z*z/q2);
	double a=sqrt(r0*r0+m*m);
	return mass/FPI/(q*pow((a+r0)*a,2)*a)*r0*(r0+2*a);
}

double rhoNFW(double R,double z){
	double m=sqrt(R*R+z*z/q2);
	return mass/(m/r0*pow(r0+m,2))/FPI;
}

/*
 *  For use as External Potential
 */
double aNFW=25.,GM_NFW=.01;
double rhoNFWext(double R,double z){
	double m=sqrt(R*R+z*z/q2);
	return GM_NFW/(m/aNFW*pow(aNFW+m,1))/FPI; // + 3.*.05/(4.*PI)*pow(1.+m*m, -2.5);
}



