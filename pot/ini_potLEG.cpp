/*
 * ini_potLEG.cpp
 *
 *  Created on: 05/mar/2014
 *      Author: morpheus
 */
#include <math.h>
#include "press.h"
#include "leg_pot3.h"

double **phil_ini, **Pr_ini, **Pr2_ini;

void intpo2_ini(double r,double *phip,int npoly){//interpolates phil to r
	if(r>=ar[nr-1]){
		for(int k=0; k<npoly; k++)
			phip[k]=phil_ini[nr-1][k]*pow(ar[nr-1]/r,2*k+1);
		} else {
		int top,bot;
		topbottom(ar,nr,r,&bot,&top);
		double db=r-ar[bot], f1=db/(ar[top]-ar[bot]);
		for(int k=0; k<npoly; k++){// linear interpolation
			phip[k]=f1*phil_ini[top][k]+(1-f1)*phil_ini[bot][k];
		}
		if(top<10){// add quadratic contributions
			int thr;
			if(f1<0.5 && bot>0) thr=bot-1; else thr=top+1;
			double dt=r-ar[top];
			double f2=dt*db/((ar[thr]-ar[top])*(ar[thr]-ar[bot])),
			f3=(ar[thr]-ar[bot])/(ar[top]-ar[bot]);
			for(int k=0; k<npoly; k++){
				phip[k]+=f2*(phil_ini[thr][k]-phil_ini[bot][k]-f3*(phil_ini[top][k]-phil_ini[bot][k]));
			}
		}
		//if (r>1e-6) printf(">> r=%f phip[0]=%f phip[1]=%f phil[0]=%f phil[1]=%f\n",r,phip[0],phip[1],phil[bot][0],phil[bot][1]);
	}
}

void intpo2_ini(double r,double *phip,double *dphip,double *d2phip,int npoly){//interpolates phil_ini to r
	if(r>=ar[nr-1]){//r larger than end of grid
		for(int k=0; k<npoly; k++){
			phip[k]=phil_ini[nr-1][k]*pow(ar[nr-1]/r,2*k+1);
			dphip[k]=-(2*k+1)*phil_ini[nr-1][k]*pow(ar[nr-1]/r,2*k+1)/r;
			d2phip[k]=(2*k+2)*(2*k+1)*phil_ini[nr-1][k]*pow(ar[nr-1]/r,2*k+1)/r/r;
		}
	} else {
		int top,bot;
		topbottom(ar,nr,r,&bot,&top);
		double db=r-ar[bot], f1=db/(ar[top]-ar[bot]);
		for(int k=0; k<npoly; k++){// linear interpolation
			phip[k]=f1*phil_ini[top][k]+(1-f1)*phil_ini[bot][k];
			dphip[k]=f1*Pr_ini[top][k]+(1-f1)*Pr_ini[bot][k];
			d2phip[k]=f1*Pr2_ini[top][k]+(1-f1)*Pr2_ini[bot][k];
		}
		if(top<10){// add quadratic contributions
			int thr;
			if(f1<0.5 && bot>0) thr=bot-1; else thr=top+1;
			double dt=r-ar[top];
			double f2=dt*db/((ar[thr]-ar[top])*(ar[thr]-ar[bot])),
			f3=(ar[thr]-ar[bot])/(ar[top]-ar[bot]);
			for(int k=0; k<npoly; k++){
				phip[k]+=f2*(phil_ini[thr][k]-phil_ini[bot][k]-f3*(phil_ini[top][k]-phil_ini[bot][k]));
				dphip[k]+=f2*(Pr_ini[thr][k]-Pr_ini[bot][k]-f3*(Pr_ini[top][k]-Pr_ini[bot][k]));
				d2phip[k]+=f2*(Pr2_ini[thr][k]-Pr2_ini[bot][k]-f3*(Pr2_ini[top][k]-Pr2_ini[bot][k]));
			}
		}
	}
}

double dPhir_ini(double r,double theta,double *d2Phidr2){
	double pol[npoly],phip[npoly],dphip[npoly],d2phip[npoly];
	double grad1=0,grad2=0;
	evenlegend(pol,cos(theta),npoly);
	intpo2_ini(r,phip,dphip,d2phip,npoly);
	for(int np=0; np<npoly; np++){
		grad1+=dphip[np]*pol[np];
		grad2+=d2phip[np]*pol[np];
	}
	(*d2Phidr2)=grad2; return grad1;
}
double dPhitheta_ini(double r, double c) {
    //returns the first derivative of the whole potential wrt spherical theta
	double dphitheta = 0;
	double dpol[npoly], phinterp[npoly];
	intpo2_ini(r, phinterp, npoly);
	for(int k=0; k<npoly; k++) dpol[k] = dlegend(c, 2*k); //fills the array dpol with derivatives of the even legendre polys
	for(int j=0; j<npoly; j++) dphitheta += phinterp[j]*dpol[j]; //this is the derivative wrt cos(theta), need to then multiply by -sin(theta)
	dphitheta *= -sqrt(1 - c*c); //im taking the positive root here - any issues with this?
	return dphitheta;
}
double d2Phitheta_ini(double r, double c) {
    //returns the second derivative of the whole potential wrt spherical theta
	double d2phitheta = 0;
	double d2pol[npoly], phinterp[npoly];
	intpo2_ini(r, phinterp, npoly);
	for(int k=0; k<npoly; k++) d2pol[k] = d2legend(c, 2*k);
	for(int j=0; j<npoly; j++) d2phitheta += phinterp[j]*d2pol[j]; //2nd derivative wrt cos(theta)
	d2phitheta *= (1 - c*c);
	d2phitheta += ( (dPhitheta(r, c)*c) / sqrt(1 - c*c) ); //add the first derivative term from the product rule and ensure the correct factor
	return d2phitheta;
}

double d2PhiR_ini(double R) {
    //returns the second derivative of Phi wrt cylindrical R in the plane
	double r = R, theta=asin(1), d2Phidr2;
	dPhir_ini(r, theta, &d2Phidr2);
	return d2Phidr2;
}
double d2Phiz_ini(double R) {
    //returns the second derivative of Phi wrt cylindrical z
	double r = R, theta=asin(1), d2;
	return (1./r)*dPhir_ini(r, theta, &d2) + (1./(r*r))*d2Phitheta_ini(r, 0);
}
double dPhiR_ini(double R, double z) {//derivative of Phi with respect to cylindrical R
	double r = sqrt(z*z + R*R), s=R/r, c=z/r, theta=acos(c), d2;
	return dPhitheta_ini(r, c)*(c/r) + dPhir_ini(r, theta, &d2)*s; //chain rule
}
double dPhiz_ini(double R, double z) {//derivative of Phi with respect to cylindrical z
	double r = sqrt(z*z + R*R), s=R/r, c=z/r, theta=acos(c), d2;
	return dPhitheta_ini(r, c)*(-s/r) + dPhir_ini(r, theta, &d2)*c; //chain rule
}

double dPhi_ini(double *x2,double *f2){
	double R=x2[0],z=x2[1];
	f2[0]=dPhiR_ini(R,z);
	f2[1]=dPhiz_ini(R,z);
	return phileg(R,z);
}

void getfreqs_ini(double R,double *kappa,double *nu,double *Omega){
	(*kappa)=sqrt(d2PhiR_ini(R)+3/R*dPhiR_ini(R,0));
	(*nu)=sqrt(d2Phiz_ini(R));
	(*Omega)=sqrt(dPhiR_ini(R,0)/R);
}
