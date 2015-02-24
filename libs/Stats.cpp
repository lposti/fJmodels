/*
 * Stats.cpp
 *
 *  Created on: Feb 24, 2015
 *      Author: morpheus
 */

#include "Grid.h"
#include "UtilsLeg.h"

// WORK HEERE PLEASE!!!
/*
void vir2(double *ci,double *wi){
	float pot=0,KRR=0,Kzz=0,WRR=0,Wzz=0,pi=acos(-1);
	for(int i=1;i<NR;i++){
		float r=.5*(ar[i-1]+ar[i]),dr=ar[i]-ar[i-1];
		for(int j=0; j<NGAUSS; j++){
			float z=r*ci[j],R=r*sqrt(1-ci[j]*ci[j]);
			double Vrot,Sigu,Sigp,Sigv,SigRz;
			float dens=ev_dens(R,z,&Vrot,&Sigu,&Sigp,&Sigv,&SigRz);
			pot+=.5*wi[j]*dens*Phi(R,z)*r*r*dr;
			KRR+=.5*wi[j]*dens*(pow(r*Sigu,2)+pow(r*Sigp,2))*dr;
			Kzz+=.5*wi[j]*dens*pow(r*Sigv,2)*dr;
			WRR-=wi[j]*dens*R*dPhiR(R,z)*r*r*dr;
			Wzz-=wi[j]*dens*z*dPhiz(R,z)*r*r*dr;
		}
	}
	KRR*=4*pi; Kzz*=4*pi; WRR*=4*pi; Wzz*=4*pi; pot*=4*pi;

	double sum=0.;
	for (int i=0;i<nr;i++) sum+=2*TPI*ar[i]*ar[i]*ev_dens(ar[i],0);  //rho_iso(ar[i]);
	printf("Mass(%3.0f): %f %f\n",ar[nr-1],pow(ar[nr-1],2)*Pr[nr-1][0],sum);
	printf("KE, PE, W/K = %f %f %f\n",KRR+Kzz,pot,pot/(KRR+Kzz));
	printf("Kxx, Wxx, Wxx/Kxx = %f %f %f\n",KRR,WRR,WRR/KRR);
	printf("Kzz, Wzz, Wzz/Kzz = %f %f %f\n",Kzz,Wzz,Wzz/Kzz);
}
*/



