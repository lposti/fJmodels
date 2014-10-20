/*
 * LM_fit.h
 *
 *  Created on: Mar 18, 2013
 *      Author: Fermani
 */

#ifndef LM_FIT_H_
#define LM_FIT_H_

#include "Structure.h"

class LM_fit{
public:
	LM_fit(double energy,int list2fit_length,double **list2fit);
	~LM_fit();
	void reset_energy(double energy);
	double surf_coeff(int ith);
	int model_f(const gsl_vector *x, void * data, gsl_vector *f);
	int model_df (const gsl_vector * x, void *data, gsl_matrix * J);
	int model_fdf (const gsl_vector * x, void *data, gsl_vector *f, gsl_matrix * J);
	void fit_qsurface_2list(int ch_l);
	int status();
	int ncoeff();
private:
	int list2fit_length_;
	double **list2fit_;
	double energy_;
	int ncoeff_;
	double fit_coeff_[9];
	int status_; // 0 for success
};


#endif /* LM_FIT_H_ */
