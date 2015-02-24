/*
 * Utils.h
 *
 *  Created on: Feb 13, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include <math.h>
#include "Grid.h"
#include "Potential.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>


#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define EPSABSROOTS 1e-8
#define EPSRELROOTS 1e-6
#define EPSABSODE   1e-10
#define EPSRELODE   1e-9
#define EPSABSINTEG 1e-5
#define EPSRELINTEG 1e-3
#define WRKSP 1000
#define PI  3.141592653589793
#define PIH 1.570796326794897
#define TPI 6.283185307179586
#define FPI 12.56637061435917
#define TINY 1e-12

/************************************************************************************
 *  Templates for arrays and matrices allocations
 ************************************************************************************/

template <typename T> inline T* arr (int const& n){
	T *t = new T[n];
	return t;
}

template <typename T> inline T** mat (int const& nx, int const& ny){
	T **t = new T*[nx];
	for (int i=0; i<nx; ++i) t[i] = new T[ny];
	return t;
}

/*************************************************************************************
 *  Templates for arrays and matrices de-allocations
 *************************************************************************************/

template <typename T> inline void delArr (T* t){
	delete [] t;
}

template <typename T> inline void delMat (T* t, int const& ny){
	for(int i=0;i<ny;i++) delete[] t[i];
	delete [] t;
}


/************************************************************************************
 *  Grid definitions
 ************************************************************************************/

// hyperbolic sine Grid
template <typename T> inline void SetGrid (T const& a_scale, T const& rmax){
	T dpsi=asinh(rmax/a_scale)/(T)(NR-1);
	for(int n=0; n<NR; n++)	ar[n]=a_scale*sinh(n*dpsi);
	ar[0]+=.25*ar[1];
}

// logarithmic Grid
template <typename T> inline void SetGrid (T const& rmax){
	T rmin=1e-4*rmax;  // minimum radius

	for (int i=0;i<NR;i++)
		ar[i] = pow(10.,log10(rmin)+(log10(rmax)-log10(rmin))*i/(NR-1));
}

template <typename T> inline T* giveLogGrid (T const& rmin, T const& rmax, int N){
	T* t = arr<T>(N);

	for (int i=0;i<N;i++)
			t[i] = pow(10.,log10(rmin)+(log10(rmax)-log10(rmin))*i/(N-1));
	return t;
}


/**********************************************************************************
 * Top-bottom and interpolation routines
 **********************************************************************************/
template <typename T> int topbottom(T *xm, int n, T x, int *botx, int *topx){
	int m, top=n-1, bot=0;
	if((xm[top]-x)*(x-xm[bot])<0.){
		if(fabs(xm[top]-x)<1e-5*fabs(x-xm[bot])){
			*topx=top; *botx=top-1; return 1;
		}else if(fabs(x-xm[bot])<1e-5*fabs(x-xm[top])){
			*topx=bot+1; *botx=bot; return 1;
		}
		else return 0;
	}
	while(top-bot>1){
		m=(top+bot)/2;
		if((xm[top]-x)*(x-xm[m])>=0) bot=m;
		else top=m;
	}
	*topx=top; *botx=bot; return 1;
}

// linear interpolation
template <typename T> T linterp(T *xm, T *ym, int np,T x){
	int top=0,bot=0;
	if(x<=xm[0]) return ym[0];
	if(x>=xm[np-1]) return ym[np-1];
	topbottom<T>(xm,np,x,&bot,&top);
	return ym[bot]+(x-xm[bot])/(xm[top]-xm[bot])*(ym[top]-ym[bot]);
}

/*************************************************************************************
 *  My-Wrapper for GSL roots
 *************************************************************************************/

template <typename T> T wrp_GSLroots(gsl_function F, T low, T high){

	const gsl_root_fsolver_type *Type;
	gsl_root_fsolver *s;

	int status=0,iter = 0;
	T Rout=0;

	Type = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (Type);
	gsl_root_fsolver_set (s, &F, low,high);

	while (status == GSL_CONTINUE || high-low > EPSABSROOTS) {
		iter++;
		status = gsl_root_fsolver_iterate (s);
		Rout  = gsl_root_fsolver_root (s);
		low    = gsl_root_fsolver_x_lower (s);
		high   = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (low, high,
				EPSABSROOTS, EPSRELROOTS);
	}

	gsl_root_fsolver_free (s);
	return Rout;
}

/*************************************************************************************
 *  My-Wrapper for GSL Integration
 *************************************************************************************/
/* 1D QAGS integration */
template <typename T> T wrp_GSLinteg(gsl_function F, T a, T b){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (WRKSP);

	T result, error;

	gsl_integration_qags (&F, a, b, EPSABSINTEG, EPSRELINTEG, WRKSP,
	                        w, &result, &error);

	gsl_integration_workspace_free (w);

	return result;
}

/*************************************************************************************
 * Find Radius if circular Orbit given Lz and E
 *************************************************************************************/
struct pELz {Potential *p; double E,Lzsq;};
template <typename T> T Rcfn(T R, void * params){//what vanishes at Rc for given Lz
	struct pELz * par = (struct pELz *)params;
	return par->Lzsq/pow(R,3)-par->p->dR(R,0);
}
template <typename T> T GetRc(T Lzsq,T R, Potential *p){
	struct pELz par = {p,0.,Lzsq};
	T Ri,Ro,dfRi,dfRo,dfRt=Rcfn<T>(R,&par);
	if(dfRt>0){//inside
		while(dfRt>0){
			Ri=R; dfRi=dfRt; R*=1.2; dfRt=Rcfn<T>(R,&par);
		} Ro=R; dfRo=dfRt;
	}else{//outside
		while(dfRt<0){
			Ro=R; dfRo=dfRt; R*=.8; dfRt=Rcfn<T>(R,&par);
		} Ri=R; dfRi=dfRt;
	}

	gsl_function F;
	F.function = &Rcfn<T>;
	F.params   = &par;

	return wrp_GSLroots(F,Ri,Ro);
}

template <typename T> T RcEfn(double R, void * params){//what vanishes at Rc for given E
	struct pELz * par = (struct pELz *)params;
	return .5*R*par->p->dR(R,0)+(*(par->p))(R,0)-par->E;
}
template <typename T> T RcE(double E, double Rc, Potential *p){
	struct pELz par = {p,E,0.};
	T R0=.8*Rc,f0=RcEfn<T>(R0,&par);
	while(f0>0){
		R0*=.8; f0=RcEfn<T>(R0,&par);
	}
	T R1=1.2*Rc, f1=RcEfn<T>(R1,&par);
	while(f1<0){
		R1*=1.2;  f1=RcEfn<T>(R1,&par);
	}
	gsl_function F;
	F.function = &RcEfn<T>;
	F.params   = &par;

	return wrp_GSLroots(F,R0,R1);
}

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

#endif /* INCLUDE_UTILS_H_ */
