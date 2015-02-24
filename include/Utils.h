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


#endif /* INCLUDE_UTILS_H_ */
