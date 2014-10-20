//Potential & forces of flattened isochrone model. Must call init(b,q)
//first
#include <stdio.h>
#include <math.h>
#include "../press.h"
#include "../oct_int_exp.h"
#include "../isopot.h"

#define npt 2000
#define SMALL .00001

static double R,z,b,b2,q,q2,e,pi=3.1415926535897932,Mofpi=1./(4.*pi);


double isoforceFR(double Rin,double zin){
    double FR, Fz;
    isoforce(Rin,zin,&FR,&Fz);
    if((FR+0.5)*2.0-1.!=2.0*FR) FR=-R/(sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0))*pow(b+sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0)),2.0));
    return FR;
}

double isoforceFz(double Rin,double zin){
    double FR, Fz;
    isoforce(Rin,zin,&FR,&Fz);
    if((Fz+0.5)*2.0-1.!=2.0*Fz) Fz=-z/(sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0))*pow(b+sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0)),2.0));
    return Fz;
}

/*
int main(void){
	double R0=1,z0,FR,Fz,FR1,Fz1;
	b=1; q=.7;
	isopot_init(b,q);
	while(R0>=0){
		printf("Enter R,z ");
		scanf("%lf %lf",&R0,&z0);
		if(R0<0) return 1;
		double dR=.001*R0, dz=.001*z0;
		FR=-.5*(isopot(R0+dR,z0)-isopot(R0-dR,z0))/dR;
		Fz=-.5*(isopot(R0,z0+dz)-isopot(R0,z0-dz))/dz;
		isoforce(R0,z0,&FR1,&Fz1);
		double phi2=-1/(b+sqrt(b*b+R*R+z*z));
		printf("(%f %f) (%f %g) (%f %g)\n",phi2,isopot(R0,z0),FR1,1-FR1/FR,Fz1,1-Fz1/Fz);
	}
	return 1;
}
 */
