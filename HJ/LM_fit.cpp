/*
 * LM_fit.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: Fermani
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "LM_fit.h"
#include "Structure.h"

int nc;
double en;

LM_fit::LM_fit(double energy,int list2fit_length,double **list2fit){
	list2fit_length_=list2fit_length;
	list2fit_=list2fit;
	if(planeorquadr==1) ncoeff_=9;
	else ncoeff_=3;
	nc=ncoeff_;
	status_=0;
	energy_=energy;
	en=energy;
}

LM_fit::~LM_fit(){}

void LM_fit::reset_energy(double energy){ energy_=energy; en=energy;}

double LM_fit::surf_coeff(int ith){ if(ith<ncoeff_) return fit_coeff_[ith]; else return 0.0;}

// Implementation of GSL multidimensional nonlinear least square fitting:
// http://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html

int model_f(const gsl_vector *x, void * data, gsl_vector *f){
	size_t n= ((struct data *)data)->n;
	double *jr = ((struct data *)data)->jr;
	double *jz = ((struct data *)data)->jz;
	double *jphi = ((struct data *)data)->jphi;
	double *sigma = ((struct data *)data)->sigma;

	double *a;
	a = new double[nc];
	for(int h=0;h<nc;h++) a[h] = gsl_vector_get(x,h);

	size_t i;
	double Yi;
	vec3 Yi_f;
	for (i = 0; i< n; i++){
		if(planeorquadr==1) Yi=a[0]*jr[i]+a[1]*jz[i]+a[2]*jphi[i]+a[3]*pow(jr[i],2.0)+a[4]*pow(jz[i],2.0)+a[5]*pow(jphi[i],2.0)+a[6]*jr[i]*jz[i]+a[7]*jz[i]*jphi[i]+a[8]*jr[i]*jphi[i];
		else Yi=a[0]*jr[i]+a[1]*jz[i]+a[2]*jphi[i];

		gsl_vector_set(f, i, (Yi-en)/sigma[i]);
	}
	return GSL_SUCCESS;
}

int model_df (const gsl_vector * x, void *data, gsl_matrix * J){
	size_t n = ((struct data *)data)->n;
	double *jr = ((struct data *)data)->jr;
	double *jz = ((struct data *)data)->jz;
	double *jphi = ((struct data *)data)->jphi;
	double *sigma = ((struct data *)data)->sigma;

	double *a;
	a = new double[nc];
	for(int h=0;h<nc;h++) a[h] = gsl_vector_get(x,h);


	size_t i;

	if(planeorquadr==1){
		for (i = 0; i < n ; i++){
			gsl_matrix_set (J, i, 0, jr[i]/sigma[i]);
			gsl_matrix_set (J, i, 1, jz[i]/sigma[i]);
			gsl_matrix_set (J, i, 2, jphi[i]/sigma[i]);
			gsl_matrix_set (J, i, 3, pow(jr[i],2.0)/sigma[i]);
			gsl_matrix_set (J, i, 4, pow(jz[i],2.0)/sigma[i]);
			gsl_matrix_set (J, i, 5, pow(jphi[i],2.0)/sigma[i]);
			gsl_matrix_set (J, i, 6, jr[i]*jz[i]/sigma[i]);
			gsl_matrix_set (J, i, 7, jphi[i]*jz[i]/sigma[i]);
			gsl_matrix_set (J, i, 8, jr[i]*jphi[i]/sigma[i]);
		}
	}
	else{
		for (i = 0; i < n ; i++){
			gsl_matrix_set (J, i, 0, jr[i]/sigma[i]);
			gsl_matrix_set (J, i, 1, jz[i]/sigma[i]);
			gsl_matrix_set (J, i, 2, jphi[i]/sigma[i]);
		}
	}
	return GSL_SUCCESS;
}

int model_fdf (const gsl_vector * x, void *data, gsl_vector *f, gsl_matrix * J){
	int model_f(const gsl_vector *x, void * data, gsl_vector *f);
	int model_df (const gsl_vector * x, void *data, gsl_matrix * J);
	model_f(x, data, f);
	model_df(x,data,J);

	return GSL_SUCCESS;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s){
	printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
			"|f(x)| = %g\n",
			iter,
			gsl_vector_get (s->x, 0),
			gsl_vector_get (s->x, 1),
			gsl_vector_get (s->x, 2),
			gsl_blas_dnrm2 (s->f));
}


void LM_fit::fit_qsurface_2list(int ch_l){
	void print_state (size_t iter, gsl_multifit_fdfsolver * s);
	int model_fdf (const gsl_vector * x, void *data, gsl_vector *f, gsl_matrix * J);
	int model_f(const gsl_vector *x, void * data, gsl_vector *f);
	int model_df (const gsl_vector * x, void *data, gsl_matrix * J);
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	int i, iter = 0;
	const size_t n = list2fit_length_;

	gsl_matrix *covar = gsl_matrix_alloc (ncoeff_, ncoeff_);
	double jr[list2fit_length_], jz[list2fit_length_], jphi[list2fit_length_], sigma[list2fit_length_];
	struct data d = {n, jr, jz, jphi, sigma};
	gsl_multifit_function_fdf f;
	double * x_init;
	x_init = new double[ncoeff_];
	for(i=0; i < ncoeff_; i++) x_init[i]=1.0;
	/*
	x_init[0]=1.0; x_init[1]=1.0; x_init[2]=1.0;
	x_init[3]=0.0; x_init[4]=0.0; x_init[5]=0.0;
	x_init[6]=1.0; x_init[7]=1.0; x_init[8]=1.0;
	*/
	gsl_vector_view x = gsl_vector_view_array (x_init, ncoeff_);
	const gsl_rng_type * type;
	gsl_rng * r;

	gsl_rng_env_setup();

	type = gsl_rng_default;
	r = gsl_rng_alloc (type);
	f.f      = &(model_f);
	f.df     = &(model_df);
	f.fdf    = &(model_fdf);
	f.n = list2fit_length_;
	f.p      = ncoeff_;
	f.params = &d;

	for (i = 0; i < list2fit_length_; i++){
		jr[i]=list2fit_[i][0];
		jz[i]=list2fit_[i][1];
		jphi[i]=list2fit_[i][2];
		sigma[i] = 0.1;
		//printf ("data: %u %g %g %g %g \n", i, jr[i], jz[i], jphi[i],sigma[i]);
	};

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, f.n, ncoeff_);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);

	print_state (iter, s);

	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		printf ("status = %s\n", gsl_strerror (status));

		print_state (iter, s);

		if (status) break;

		status = gsl_multifit_test_delta (s->dx, s->x,1e-8, 1e-8);
	}
	while (status == GSL_CONTINUE && iter < ch_l);

	gsl_multifit_covar (s->J, 0.0, covar);

	 printf ("status = %s\n", gsl_strerror (status));
	 status_=status; // 0 for success
	 if(ncoeff_>3){
		 for ( i = 0; i < 9; i++){
			 fit_coeff_[i]=gsl_vector_get(s->x, i);
		 }
	 }
	 else{
		 for ( i = 0; i < ncoeff_; i++){
			 fit_coeff_[i]=gsl_vector_get(s->x, i);
		 }
		 for (i = ncoeff_; i< 9; i++){
			 fit_coeff_[i]=0.0;
		 }
	 }

	 gsl_multifit_fdfsolver_free (s);
	 gsl_matrix_free (covar);
	 gsl_rng_free (r);
}

int LM_fit::status(){ return status_;}

int LM_fit::ncoeff(){ return ncoeff_;}
