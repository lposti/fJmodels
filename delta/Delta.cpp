/*
 * Delta.cpp
 *
 *  Created on: Feb 16, 2015
 *      Author: morpheus
 */

#include "Delta.h"
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "Utils.h"


/* Constructor ********************************************************/
Delta::Delta(struct eqsPar * par_in) {
	par = par_in;
	v0=0.; cv=0.; sv=0.; R0=0.;
	Ri=0.; zi=0.; sqDR=0.; D2=0;
	size=200;

	Rc=RcE<double>(par->E,par->R0,&(par->p));
}

/* v^2 */
inline double vsq(struct eqsPar * par,double x){
	return 2*(par->E-(par->p(x,0.)+.5*pow(par->Lz/x,2)));
}

/* Equations of Motion in Cylindrical coords. */
int eqsMotion(double t,const double y[], double f[],void *parameters){
        struct eqsPar *par = (struct eqsPar *)(parameters);

        f[0] = y[2];
        f[1] = y[3];
        f[2] = -par->p.dR(y[0],y[1])+pow(par->Lz,2.0)/pow(y[0],3.0);
        f[3] = -par->p.dz(y[0],y[1]);

        return GSL_SUCCESS;
}

/*
 *  Computes the orbit from the equatorial plane with pR=0 and pz=vmax.
 *  Uses GSL's RK8PD integrator with an adaptive step-finding controller
 *  that keeps fixed the accuracy of the integration.
 *  Stores (R,z,pR,pz) in constant sized arrays (size). The exit condition
 *  must be passed as a (condInterf*) functor.
 *
 *  Writes the integrated orbit in 'orb.dat'
 */
unsigned Delta::evolveOrbit(double x0, double * R, double * z, double * pR, double * pz,
							const unsigned size, condInterf * cond){

	double steph=2e-3, y[4], t=0;
	y[0]= x0; y[1]= 0.; y[2]=0.; y[3]= sqrt(vsq(par,x0));

	const gsl_odeiv2_step_type * T
	    = gsl_odeiv2_step_rk8pd;

	  gsl_odeiv2_step * s
	    = gsl_odeiv2_step_alloc (T, 4);
	  gsl_odeiv2_control * c
	    = gsl_odeiv2_control_y_new (EPSABSODE, EPSRELODE);
	  gsl_odeiv2_evolve * e
	    = gsl_odeiv2_evolve_alloc (4);

	  gsl_odeiv2_system sys = {eqsMotion, NULL, 4, par};

    FILE* fp=fopen("orb.dat","w");
    unsigned i=0;
    R[i]=y[0]; z[i]=y[1]; pR[i]=y[2]; pz[i]=y[3];

    while ((*cond)(y,i,size)){
    	gsl_odeiv2_evolve_apply (e, c, s,&sys,
    	                         &t, 5*x0/sqrt(vsq(par,x0)), // 5 times the typical timescale
  	                             &steph, y);                 // is a typical upper-limit

    	 i++; R[i]=y[0]; z[i]=y[1]; pR[i]=y[2]; pz[i]=y[3]; // update variables

    	 fprintf(fp,"%f %f %f\n",R[i],z[i],steph);
    }
    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
    fclose(fp);

    return i;
}

/* Get the difference btw. R0 and crossing radius */
double Delta::getDR(double x0){

	double *R = new double [size],*z = new double [size],*pR = new double [size],*pz = new double [size];
	/* Explicit casting from functor class to interface */
	Cond_zgt0 * cz = new Cond_zgt0;
	condInterf * cz_b = (condInterf *)cz;
	int N = evolveOrbit(x0,R,z,pR,pz,size,cz_b);

	if (fabs(z[N-1]-z[N])>0 and fabs(R[N-1]-R[N]))
		return x0 - (R[N]-z[N]*(R[N-1]-R[N])/(z[N-1]-z[N]));
	else return x0 - R[N];
}

/* Wrapper for getDR (used by GSL functions) */
double wrp_getDR(double x0, void * params){
	Delta * dd = (Delta *)params;
	return dd->getDR(x0);
}

/* Derivative of effective potential */
inline double Delta::dPhiEff(const double x){ return par->p.dR(x,0)-pow(par->Lz,2)/pow(x,3);}

/* Returns the nearest x0 for which v^2>0 */
inline double Delta::ensureVsq(double x0){
	int count=0;
	while(vsq(par,x0)<0){
		if (dPhiEff(x0)<0) x0*=1.1;
		else x0*=.9;
		count++; if(count>50){ printf("vsq: %f %f %f %f %f\n",x0,par->E,par->Lz,par->p.dR(x0,0),vsq(par,x0)); exit(0);}
	}
	return x0;
}
/*
 *  Find R0 of Jr=0 orbit (at given E).
 *  First bracket the root, then uses GSL root finder (Brent)
 */
double Delta::getJr0(){
	double x0=Rc,x1;
	double dx0,dx1;
	double A=.8;
	int k=0;

	x0=ensureVsq(x0);   dx0=getDR(x0);
	x1=ensureVsq(x0/A); dx1=getDR(x1);
	while(dx0*dx1>0){//we haven't bracketed the root
		if(fabs(dx1)<fabs(dx0)){//more of same
			x0=x1; dx0=dx1;
		}else{// back off
			A=pow(A,-0.9);
		}
		x1=ensureVsq(x0*A); dx1=getDR(x1);
		k++;
		if(k==21){
			printf("get_closed: %f %f %g %g %g\n",x0,x1,dx0,dx1,A);
			if(fabs(dx0)<1.e-7) return -x0;
			else{
				printf("problem in get_closed(): %f %f %f %f\n",x0,x1,dx0,dx1);
				return -1;
			}
		}
	}

	double low=MIN(x0,x1),high=MAX(x0,x1);

	gsl_function F;

	F.function = &wrp_getDR;
	F.params   = this;

	if (!(low<high)) printf("Delta::getJr0: %f %f\n",low,high);
	return wrp_GSLroots<double>(F,low,high);
}

/*
 *  Derivatives of s^2 w.r.t. v and Delta2
 */

double Delta::ds2dv(const double v){
	return 2*(zi*sqDR*sin(v)-cos(v)*(Ri*R0+D2*sin(v)));
}

/*
 *  Wrapper needed for GSL - C++ classes integration:
 *  in the void * params of the gsl_function pass also
 *  a pointer to the class (use this in the params list).
 */
double wrp_ds2dv(double v, void * params){
	Delta *dd = (Delta *) params;
	return dd->ds2dv(v);
}

/*
 *  Get angle at which the closest point of the current
 *  ellipse is: given D2, finds v which minimizes deviation
 *  from the given point (Ri,zi)
 */
double Delta::getv(){

    double low=0,high=PIH;

    if      (fabs(ds2dv(low)) <EPSABSROOTS) v0=low;
    else if (fabs(ds2dv(high))<EPSABSROOTS) v0=high;
		else {

			gsl_function F;

			F.function = &wrp_ds2dv;
			F.params   = this;

			if (!(low<high)) printf("Delta::getv: %f %f\n",low,high);
			v0 = wrp_GSLroots(F,low,high);
		}

	sv=sin(v0); cv=cos(v0);
	return v0;
}

double Delta::ds2dD2TOT(double D2h, void * params){

	struct Rzptr * p_Rz = (struct Rzptr *) params;
	double *Rpoint = p_Rz->R, *zpoint = p_Rz->z;
	int N = p_Rz->N;

	D2 = D2h; sqDR=sqrt(D2+R0*R0);

	double sum=0.;
	for (int i=0; i<N ;i++){
		/* Given a point on the orbit, find the v of the point
		 * in the ellipse (with a given Delta^2) closest to
		 * that point.
		 */
		Ri = Rpoint[i]; zi = zpoint[i];
		getv();

		/* Then compute the sum of the least squares of those
		 * points (in the same ellipse at fixed Delta^2)
		 */
		double dvdD2 = .5*sv*(zi/sqDR-2*cv)/(cv*zi*sqDR+sv*Ri*R0-D2*(cv*cv-sv*sv));
		sum+=(cv-zi/sqDR)*cv+ds2dv(v0)*dvdD2;
	}
	return sum;
}

/* Wrapper for GSL functions */
double wrp_ds2dD2TOT(double D2h, void * params){
	struct Rzptr * par = (struct Rzptr *) params;
	return par->dd->ds2dD2TOT(D2h,params);
}

/*
 *  Get best focal distance of ellipsoidal coords.
 *  First finds the orbit with Jr=0 at that energy, then
 *  finds the best ellipse which minimizes:
 *  (R-R0*sin(v))^2 + (z-sqrt(R0^2+D^2)*cos(v))^2
 */
double Delta::getBestDelta(){

	R0=getJr0();

	double *R = new double [size],*z = new double [size],*pR = new double [size],*pz = new double [size];

	/* Explicit casting from functor class to interface */
	Cond_pzgt0 * cpz = new Cond_pzgt0;
	condInterf * cpz_b = (condInterf *)cpz;
	int N = evolveOrbit(R0,R,z,pR,pz,size,cpz_b);

	struct Rzptr p_Rz = {R,z,N,this};

	double low=-0.001,high=5, D2out=0;
	double Dmin=ds2dD2TOT(low,&p_Rz), Dmax=ds2dD2TOT(high,&p_Rz);

	while(Dmin*Dmax>0 and fabs(Dmax)>1e-4 and fabs(Dmin)>1e-4){
		if(Dmax<0){
			high*=5; Dmax=ds2dD2TOT(high,&p_Rz);
		}else{
			low-=.1; Dmin=ds2dD2TOT(low,&p_Rz);
		}
	}

	if      (fabs(ds2dD2TOT(low,&p_Rz))  < EPSABSROOTS) D2out = low;
	else if (fabs(ds2dD2TOT(high,&p_Rz)) < EPSABSROOTS) D2out = high;
		else {

			gsl_function F;

			F.function = &wrp_ds2dD2TOT;
			F.params   = &p_Rz;

			if (!(low<high)) printf("Delta::getBestDelta: %f %f\n",low,high);
			D2out = wrp_GSLroots(F,low,high);
		}

	delete [] R;  delete [] z;
	delete [] pR; delete [] pz;

	/*TEST*/
	return D2out;

	if (D2out<0) return -sqrt(fabs(D2out));
	else return sqrt(D2out);
	//else return (((0.01)>(sqrt(D2out))?(0.01):(sqrt(D2out)))); // MAX between 0 and sqrt(D2out)
}
