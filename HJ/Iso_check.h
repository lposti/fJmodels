/*
 * Iso_check.h
 *
 *  Created on: Mar 16, 2013
 *      Author: Fermani
 */

#ifndef ISO_CHECK_H_
#define ISO_CHECK_H_

class Iso_check{
public:
	Iso_check(int npar,double *coeff);
	~Iso_check();
	double Jr(double E,double Jz,double Jphi);
	double Jphi(double E,double Jr,double Jz);
	double E(double Jr,double Jz,double Jphi);
	double L(double E, double Jr);
	double DeltaJr(double **scoeff,int sindex,double Jz,double Jphi);
	double DeltaE(double **scoeff,int sindex,double Jz,double Jphi);
	double rho(double r);
private:
	int npar_;
	double a_[5];
};


#endif /* ISO_CHECK_H_ */
