/*
 * Stats.cpp
 *
 *  Created on: Feb 24, 2015
 *      Author: L. Posti
 */

#include "Potential.h"
#include "Grid.h"
#include "UtilsLeg.h"

void vir2(Potential *p){

	double ci[NGAUSS],wi[NGAUSS];
	gauleg<double>(0,1,ci,wi);

	double pot=0,KRR=0,Kzz=0,WRR=0,Wzz=0;
	double Vrot,Sigu,Sigp,Sigv,SigRz;

	for(int i=1;i<NR;i++){
		double r=.5*(ar[i-1]+ar[i]),dr=ar[i]-ar[i-1];
		for(int j=0; j<NGAUSS; j++){
			double sij = sqrt(1-ci[j]*ci[j]);
			double z=r*ci[j],R=r*sij;
			double dens=ev_dens(R,z,&Vrot,&Sigu,&Sigp,&Sigv,&SigRz);
			pot+=.5*wi[j]*dens*(*p)(R,z)*r*r*dr;
			KRR+=.5*wi[j]*dens*(pow(r*Sigu,2)+pow(r*Sigp,2))*dr;
			Kzz+=.5*wi[j]*dens*pow(r*Sigv,2)*dr;
			WRR-=wi[j]*dens*R*p->dR(R,z)*r*r*dr;
			Wzz-=wi[j]*dens*z*p->dz(R,z)*r*r*dr;
		}
	}
	KRR*=4*PI; Kzz*=4*PI; WRR*=4*PI; Wzz*=4*PI; pot*=4*PI;

	double sum=0.;
	for (int i=0;i<NR;i++) sum+=2*TPI*ar[i]*ar[i]*ev_dens<double>(ar[i],0);  //rho_iso(ar[i]);
	printf("Mass(%3.0f): %f %f\n",ar[NR-1],pow(ar[NR-1],2)*Pr[NR-1][0],sum);
	printf("KE, PE, W/K = %f %f %f\n",KRR+Kzz,pot,pot/(KRR+Kzz));
	printf("Kxx, Wxx, Wxx/Kxx = %f %f %f\n",KRR,WRR,WRR/KRR);
	printf("Kzz, Wzz, Wzz/Kzz = %f %f %f\n",Kzz,Wzz,Wzz/Kzz);
}




