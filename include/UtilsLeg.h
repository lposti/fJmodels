/*
 * UtilsLeg.h
 *
 *  Created on: Feb 24, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_UTILSLEG_H_
#define INCLUDE_UTILSLEG_H_

#include "Grid.h"
#include "Utils.h"

/*************************************************************************************
 * Potential stuff
 *************************************************************************************/
#define EPS 3.0e-11
template <typename T> void gauleg(T x1, T x2, T *x, T *w){
	int n=NGAUSS;
	int m=(n+1)/2;
	T xm=0.5*(x2+x1), xl=0.5*(x2-x1), pp, z1;
	for (int i=0;i<m;i++) {
		T z=cos(3.141592654*(i+0.75)/(n+0.5));
		do {
			T p1=1, p2=0, p3;
			for (int j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*(j+1)-1)*z*p2-j*p3)/(T)(j+1);
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z; z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[n-i-1]=xm+xl*z;
		w[i]=2*xl/((1-z*z)*pp*pp);
		w[n-i-1]=w[i];
	}
}
#undef EPS


template <typename T> void legend(T *allpol, T const c, int npoly) {
    // evaluates the legendre polys up to l = npoly at c ------------
	npoly++;
	allpol[0] = 1; if(npoly<2) return;
	allpol[1] = c;
	for(int i = 2; i < npoly; i++)
		allpol[i] = ((2*i-1)*c*allpol[i-1] - (i-1)*allpol[i-2]) /(double)i;
}
template <typename T> T dlegend(T const c, int n) {
    // evaluates the derivative of the legendre polynomial of order n
	if(n == 0) return 0;
	double allpol[n+1];
	legend(allpol, c, n);
	return (n*allpol[n-1] - n*c*allpol[n]) / (1 - c*c);
}


// evaluates even legendre polys up to l=2*(npoly-1) at c -------------
template <typename T> void evenlegend(T *pol,T c){
	int npoly=NPOLY;
	double c2=c*c;
	pol[0]=1; if(npoly<2) return;
	pol[1]=1.5*c2-.5;
	for(int np=2; np<npoly; np++){
		int l=2*(np-1), l2=2*l;
		pol[np]=-pol[np-2]*l*(l-1)/(T)((l2+1)*(l2-1))+
			  pol[np-1]*(c2-(l2*l+l2-1)/(T)((l2-1)*(l2+3)));
		pol[np]*=(l2+1)*(l2+3)/(T)((l+1)*(l+2));
	}
}

template <typename T> void fillPoly(T *si, T* ci, T* wi, T* pol, T **poly){
	// npoly coeffs for Gauss integration over angles ----------------------------
	int ngauss=NGAUSS,npoly=NPOLY;

	for(int i=0; i<ngauss; i++)
			si[i]=sqrt(1-ci[i]*ci[i]);

	for(int i=0; i<ngauss; i++){
		evenlegend<T>(pol,ci[i]);
		for(int np=0; np<npoly; np++) poly[np][i]=2*pol[np]*wi[i];
	}
}

// interpolates Legendre coefficients for the potential
template <typename T> void intpo2(T r, T *phip){
	int npoly=NPOLY;
	if(r>=ar[NR-1]){
		for(int k=0; k<npoly; k++)
			phip[k]=phil[NR-1][k]*pow(ar[NR-1]/r,2*k+1);

	} else {
		/*
		 * 02/09/14: added check on minimum radius, so that I won't get seg. faults.
		 * 			 to be checked its consistency anyway
		 */
		int top=0,bot=0; r=MAX(r,ar[0]);
		topbottom<T>(ar,NR,r,&bot,&top);
		T db=r-ar[bot], f1=db/(ar[top]-ar[bot]);
		for(int k=0; k<npoly; k++){// linear interpolation
			phip[k]=f1*phil[top][k]+(1-f1)*phil[bot][k];
		}
		if(top<10){// add quadratic contributions
			int thr;
			if(f1<0.5 && bot>0) thr=bot-1; else thr=top+1;
			T dt=r-ar[top];
			T f2=dt*db/((ar[thr]-ar[top])*(ar[thr]-ar[bot])),
			f3=(ar[thr]-ar[bot])/(ar[top]-ar[bot]);

			for(int k=0; k<npoly; k++)
				phip[k]+=f2*(phil[thr][k]-phil[bot][k]-f3*(phil[top][k]-phil[bot][k]));
		}
	}
}

// interpolates Legendre coefficients for the potential and its first and second derivative
template <typename T> void intpo2(T r,T *phip,T *dphip,T *d2phip){
	unsigned nr=NR,npoly=NPOLY;
	if(r>=ar[nr-1]){//r larger than end of grid
		for(int k=0; k<npoly; k++){
			phip[k]=phil[nr-1][k]*pow(ar[nr-1]/r,2*k+1);
			dphip[k]=-(2*k+1)*phil[nr-1][k]*pow(ar[nr-1]/r,2*k+1)/r;
			d2phip[k]=(2*k+2)*(2*k+1)*phil[nr-1][k]*pow(ar[nr-1]/r,2*k+1)/r/r;
		}
	} else {
		int top=0,bot=0;
		topbottom(ar,nr,r,&bot,&top);
		double db=r-ar[bot], f1=db/(ar[top]-ar[bot]);
		for(int k=0; k<npoly; k++){// linear interpolation
			phip[k]=f1*phil[top][k]+(1-f1)*phil[bot][k];
			dphip[k]=f1*Pr[top][k]+(1-f1)*Pr[bot][k];
			d2phip[k]=f1*Pr2[top][k]+(1-f1)*Pr2[bot][k];
		}
		if(top<10){// add quadratic contributions
			int thr;
			if(f1<0.5 && bot>0) thr=bot-1; else thr=top+1;
			double dt=r-ar[top];
			double f2=dt*db/((ar[thr]-ar[top])*(ar[thr]-ar[bot])),
			f3=(ar[thr]-ar[bot])/(ar[top]-ar[bot]);
			for(int k=0; k<npoly; k++){
				phip[k]+=f2*(phil[thr][k]-phil[bot][k]-f3*(phil[top][k]-phil[bot][k]));
				dphip[k]+=f2*(Pr[thr][k]-Pr[bot][k]-f3*(Pr[top][k]-Pr[bot][k]));
				d2phip[k]+=f2*(Pr2[thr][k]-Pr2[bot][k]-f3*(Pr2[top][k]-Pr2[bot][k]));
			}
		}
	}
}

/*
 * Calculates Legendre poly coeffs Ql for the quantity Q
 */
template <typename T> void get_Ylm(T **Q,T **Ql){

	int nr=NR,ngauss=NGAUSS,npoly=NPOLY;
	double ci[ngauss],wi[ngauss],pol[npoly],**poly;

	poly=mat<T>(npoly,ngauss);
	gauleg<T>(0,1,ci,wi);

	/* store Leg. polys in pol */
	for(int ng=0; ng<ngauss; ng++)
		evenlegend<T>(pol,ci[ng]);

	for(int ng=0; ng<ngauss; ng++)
		for(int np=0; np<npoly; np++)
			poly[np][ng]=2*pol[np]*wi[ng];

	/* first set to zeros */
	for(int n=1; n<nr; n++)
		for(int np=0; np<npoly; np++)
			Ql[n][np]=0;

	/* compute */
	for(int n=1; n<nr; n++)
		for(int np=0; np<npoly; np++)
			for(int ng=0; ng<ngauss; ng++)
				Ql[n][np]+=poly[np][ng]*Q[n][ng];


	Ql[0][0]=2*Q[0][0];
	for(int np=1; np<npoly; np++) Ql[0][np]=0;

	delMat(poly,npoly);
}

/*
 *  Compute Legendre poly coeffs for rho,sigR,sigp,sigz,sigRz
 */
template <typename T> void get_Ylm(T **rho,T **rhL,T **vrot,T **vrotL,
								   T **sigR,T **sigRL,T **sigp,T **sigpL,
								   T **sigz,T **sigzL,T **sigRz,T **sigRzL){

	int nr=NR,ngauss=NGAUSS,npoly=NPOLY;
	double ci[ngauss],wi[ngauss],pol[npoly],**poly;

	poly=mat<T>(npoly,ngauss);
	gauleg<T>(0,1,ci,wi);

	/* store Leg. polys in pol */
	for(int ng=0; ng<ngauss; ng++)
		evenlegend<T>(pol,ci[ng]);

	for(int ng=0; ng<ngauss; ng++)
		for(int np=0; np<npoly; np++)
			poly[np][ng]=2*pol[np]*wi[ng];

	/* first set to zeros */
	for(int n=0; n<nr; n++)
		for(int np=0; np<npoly; np++){
			rhL[n][np]=0;   vrotL[n][np]=0;
			sigRL[n][np]=0; sigpL[n][np]=0;
			sigzL[n][np]=0; sigRzL[n][np]=0;
		}

	/* compute */
	for(int n=0; n<nr; n++)
		for(int np=0; np<npoly; np++)
			for(int ng=0; ng<ngauss; ng++){
				rhL[n][np]+=poly[np][ng]*rho[n][ng];
				vrotL[n][np]+=poly[np][ng]*vrot[n][ng];
				sigRL[n][np]+=poly[np][ng]*sigR[n][ng];
				sigpL[n][np]+=poly[np][ng]*sigp[n][ng];
				sigzL[n][np]+=poly[np][ng]*sigz[n][ng];
				sigRzL[n][np]+=poly[np][ng]*sigRz[n][ng];
			}


	rhL[0][0]=2*rho[0][0]; vrotL[0][0]=2*vrot[0][0]; sigRL[0][0]=2*sigR[0][0];
	sigpL[0][0]=2*sigp[0][0]; sigzL[0][0]=2*sigz[0][0]; sigRzL[0][0]=2*sigRz[0][0];

	for(int np=1; np<npoly; np++){
		rhL[0][np]=0; vrotL[0][np]=0; sigRL[0][np]=0;
		sigpL[0][np]=0; sigzL[0][np]=0; sigRzL[0][np]=0;
	}

	delMat(poly,npoly);
}

/************************************************************************************
 *  Density stuff
 ************************************************************************************/

template <typename T> void int_dens(T r,T *densp){
	if(r>ar[NR-1]){
		for(int i=0; i<NPOLY; i++){
			densp[i]=0;
		}
	}else{
		int top,bot;
		topbottom(ar,NR,r,&bot,&top);
		double f=(r-ar[bot])/(ar[top]-ar[bot]);
		for(int i=0; i<NPOLY; i++){
			densp[i]=f*rhl[top][i]+(1-f)*rhl[bot][i];
		}
	}
}

template <typename T> void int_dens(T r,T *densp,T *Vrotp,
									T *sigup,T *sigpp,T *sigvp){
	if(r>ar[NR-1]){
		for(int i=0; i<NPOLY; i++){
			densp[i]=0; Vrotp[i]=0; sigup[i]=0; sigpp[i]=0; sigvp[i]=0;
		}
	}else{
		int top,bot;
		topbottom<T>(ar,NR,r,&bot,&top);
		float f=(r-ar[bot])/(ar[top]-ar[bot]);
		for(int i=0; i<NPOLY; i++){
			densp[i]=f*rhl[top][i]+(1-f)*rhl[bot][i];
			sigup[i]=f*sigRl[top][i]+(1-f)*sigRl[bot][i];
			sigpp[i]=f*sigpl[top][i]+(1-f)*sigpl[bot][i];
			sigvp[i]=f*sigzl[top][i]+(1-f)*sigzl[bot][i];
			//Vrotp[i]=f*vbarl[top][i]+(1-f)*vbarl[bot][i];
		}
	}
}

template <typename T> T ev_dens(T R,T z){//returns dens from rhl
	double pol[NPOLY],densp[NPOLY];
	double c=z/sqrt(R*R+z*z), r=sqrt(R*R+z*z);
	evenlegend<T>(pol,c); int_dens<T>(r,densp);
	double dens=.5*densp[0];
	for(int np=1; np<NPOLY; np++){
		double fac=.5*(4*np+1);
		dens+=fac*densp[np]*pol[np];
	}
	return dens;
}



#endif /* INCLUDE_UTILSLEG_H_ */
