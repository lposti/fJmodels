/*
 * computeJtheta.cpp
 *
 *  Created on: 27/giu/2014
 *      Author: LP
 */

#include "surf_lib.h"
#include "Types.h"
#include "Pi.h"
#include "../oct_int_exp.h"
#include "tools.h"
#include "eqs_of_motion.h"
#include "../press.h"

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

using namespace std;

#define VERBOSE


double gslInterpInteg( double * x, double * p, int n ){

	gsl_set_error_handler_off();

	const gsl_interp_type *interpType;
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	if (n>2) interpType = gsl_interp_cspline;
	else interpType = gsl_interp_linear;

	gsl_interp * interp = gsl_interp_alloc(interpType, n);
	if (gsl_interp_init(interp, x, p, n) != GSL_SUCCESS) {
		for (int i=0; i<n; i++) if (x[i]>=x[i+1]) cout << i << " " << x[i] << " " << p[i] << " " << x[i+1] << " " << p[i+1]<< endl;
		exit(1);
	}
	double res = gsl_interp_eval_integ (interp, x, p, x[1], x[n-1], acc);

	gsl_interp_free (interp);
	gsl_interp_accel_free (acc);

	return res;
}


int evolveOrbitSpherical(eveq_p * par, double rini,
							double * r, double * Pr,
							double * theta, double * Ptheta,
							double * thetaMin, double * rfin,
							int * n, int * nr, int * flag,
							double steph, int size) {

	(*flag)=0;
	int i=0,jr=0,kr=0; *n=0;
	double y[4],eps=1e-6,t=0,ti=0;
	double ptheta = -sqrt(pow(par->L,2)-pow(par->Lz,2)),
		pr = sqrt(2.*(par->energy-par->p->phi(rini,0))-pow(par->L/rini,2));

#ifdef VERBOSE
	cout << "=== check on Pr: " << ( isnan(pr)==1 ? 0 : pr ) << endl;
#endif

	y[0]= rini; y[1]= Pi/2; y[2]= ( isnan(pr)==1 ? 0 : pr ); y[3]= ptheta;
	gsl_odeiv2_system sys = {eqs_of_motion_spherical, NULL, 4, par}; // GSL integrator tested and functioning
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, steph, eps, eps);

	theta[0]=y[1]; Ptheta[0]=y[3]; r[0]=y[0]; Pr[0]=y[2];
	*thetaMin=Pi; double rmax=rini,rmin=rini;

	bool thetaNotPi2 = true, evalThetaNotPi2 = false;
	FILE* orbp = fopen("orbit.dat","w");
	/*
	 * exit conditions as follows:
	 * 1) max. number of iterations "size" reached
	 * 2) theta ~ Pi/2
	 * the second condition is evaluated only after
	 * theta has dropped below Pi/3, before that it
	 * is defaulted to false.
	 */
	while (i<size && thetaNotPi2 ){
		if (evalThetaNotPi2 or y[1]<Pi/3) {
			evalThetaNotPi2 = true;
			thetaNotPi2 = fabs(y[1]-Pi/2)>1e-2;
		}
		i++;
		gsl_odeiv2_driver_apply (d, &ti, t, y);

		/* theta is in [0,Pi] */
		if (par->Lz==0. and y[1]<1e-2) {
			/* if we are at the inversion of motion
			 * then invert sign of ptheta (only needed
			 * for polar orbits in spherical potentials
			 * since ptheta=const)
			 */
#ifdef VERBOSE
			cout << "INVERSION OF MOTION @ (theta,Ptheta)=" << y[1] << "," << y[3] << endl;
#endif
			y[3]=-y[3];
		}
		y[1]=fabs(y[1]);

		if (y[1]<*thetaMin) {*thetaMin=y[1]; (*n)++;}
		if (y[0]<rmin) {rmin=y[0]; jr++;}
		if (y[0]>rmax) {rmax=y[0]; kr++;}
		if (y[1]>.6*Pi) (*flag)=1;
		Ptheta[i]=y[3]; theta[i]=y[1]; r[i]=y[0]; Pr[i]=y[2];
		ti=t;
		t=t+steph;
		fprintf(orbp,"%e %e %e %e\n",y[0],y[1],y[2],y[3]);
	}
	fclose(orbp);

	gsl_odeiv2_driver_free (d);
	(*nr)= ( jr>0 ? jr : kr );
	(*rfin) = r[i];
	return i;
}

double getStep(eveq_p * par) {
	double maxE = par->p->phi(0.05,0.05), minE = par->p->phi(75,75);
	double maxStep = 1e-3, minStep = 5e-5;

	// linear interpolation
	return minStep + (maxStep - minStep)/(minE - maxE)*(par->energy-maxE);
}

vec2 computeJtheta(eveq_p * par, double Rc, int Jrint) {

	int size=5e4,nn,nnr,fin=size,flag=0,n=0,j=0,m=0;
	double theta[size], r[size], Ptheta[size], Pr[size], thetaMin, rfin;

	/* get proper step for orbital integration */
	double step = getStep(par);

	fin=evolveOrbitSpherical(par, Rc, r, Pr, theta, Ptheta, &thetaMin, &rfin, &nn, &nnr, &flag, step, size);
	while (fin==size){
#ifdef VERBOSE
		cout << endl << "\t @@ old/new step: " << step;
#endif
		if (flag==1) step/=2;
		else step*=2;
#ifdef VERBOSE
		cout << "/" << step << endl;
#endif
		fin=evolveOrbitSpherical(par, Rc, r, Pr, theta, Ptheta, &thetaMin, &rfin, &nn, &nnr, &flag, step, size);
	}

	cout << "--> Initial radius: " << Rc <<
			",  Final radius: " << rfin << endl;

	for (int i=0; i<nn+1; i++) if (i>0 and theta[nn-i]!=theta[nn-i-1]) n++;
	for (int i=nn+1; i<fin; i++) if (theta[fin-i]!=theta[fin-i-1]) m++;

	double Ptheta2[n],theta2[n],Ptheta2Inv[m],theta2Inv[m]; double res,res2;

	/* orbit on the direct path: from Pi/2 to thetaMin */
	Ptheta2[j] = Ptheta[nn]; theta2[j]  = theta[nn]; j++;
	for (int i=1; i<nn+1; i++){
		if (theta[nn-i]!=theta[nn-i-1]){
			Ptheta2[j] = Ptheta[nn-i];
			theta2[j]  = theta[nn-i];

			j++;
		}
	}
	res = gslInterpInteg(theta2,Ptheta2,n);


	/* orbit on the return path: from thetaMin to Pi/2 */
	j=0;
	for (int i=nn+1; i<fin; i++){
		if (theta[fin-i]!=theta[fin-i-1]){
			Ptheta2Inv[j] = Ptheta[i];
			theta2Inv[j]  = theta[i];

			j++;
		}
	}
	res2 = gslInterpInteg(theta2Inv, Ptheta2Inv, m);


	/* return value: array of two doubles */
	vec2 J;
	J[0] = -1/Pi*res+1/Pi*res2;
#ifdef VERBOSE
	cout << "Computed Jtheta: " << J[0] << ", compare with L+|Lz|: " << par->L-fabs(par->Lz) << endl;
#endif


	/*
	 * Compute Jr integral, if needed
	 */
	if (Jrint) {
		int nr=0,mr=0;
		for (int i=0; i<nnr+1; i++) if (i>0 and r[nnr-i]!=r[nnr-i-1]) nr++;
		for (int i=nnr+1; i<fin; i++) if (r[fin-i]!=r[fin-i-1]) mr++;

		double Pr2[nr],r2[nr],Pr2Inv[mr],r2Inv[mr]; double resR,res2R;


		j=0;
		/* orbit on the direct path: from thetaMin to Pi/2 */
		if (r[0]>r[1]){
			/* from end of array to beginning */
			Pr2[j] = Pr[nnr]; r2[j]  = r[nnr]; j++;
			for (int i=1; i<nnr+1; i++){
				if (r[nnr-i]!=r[nnr-i-1]){
					Pr2[j] = Pr[nnr-i];
					r2[j]  = r[nnr-i];

					j++;
				}
			}
			resR = gslInterpInteg(r2,Pr2,nr);
		} else {
			/* from beginning to end of array */
			for (int i=0; i<nnr+1; i++){
				if (r[i]!=r[i-1]){
					Pr2[j] = Pr[i];
					r2[j]  = r[i];
					j++;
				}
			}
			resR = gslInterpInteg(r2,Pr2,nr);
		} //end if direct path


		j=0;
		/* orbit on the return path: from thetaMin to Pi/2 */
		if (r[nnr+1]>r[nnr+2]){
			/* from end of array to beginning */
			Pr2[j] = Pr[fin]; r2[j]  = r[fin]; j++;
			for (int i=nnr+1; i<fin; i++){
				if (r[fin-i]!=r[fin-i-1]){
					Pr2Inv[j] = Pr[fin-i+nnr+1];
					r2Inv[j]  = r[fin-i+nnr+1];
					j++;
				}
			}
			res2R = gslInterpInteg(r2Inv,Pr2Inv,mr);
		} else {
			/* from beginning to end of array */
			for (int i=nnr+1; i<fin; i++){
				if (r[i]!=r[i-1]){
					Pr2Inv[j] = Pr[i];
					r2Inv[j]  = r[i];

					j++;
				}
			}
			res2R = gslInterpInteg(r2Inv,Pr2Inv,mr);
		} //end if return path


		J[1] = 1/Pi*(resR+res2R);
#ifdef VERBOSE
		cout << "Computed Jr    : " << J[1] << endl;
#endif

	} //end if (Jrint)

	return J;
}


/*
 * APIs: getJtheta, if only Jtheta is needed...  return value: double
 * 		 getJthetaJr, if also Jr is needed.....  return value: vec2
 */
double getJtheta (eveq_p * par, double Rc) {
	vec2 Jtheta = computeJtheta(par,Rc,0);
	return Jtheta[0];
}

vec2 getJthetaJr (eveq_p * par, double Rc) {
	return computeJtheta(par,Rc,1);
}
