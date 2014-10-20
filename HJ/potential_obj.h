/*
 * potential.h
 *
 *  Created on: Mar 15, 2013
 *      Author: Fermani
 */

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "falPot.h"

class potential_obj{
public:
	potential_obj(int pot_index,int npar,double *coeff,char *fname);
	/*
	 * 14/05/14
	 * LP edit: made the following methods virtuals so that they can be overridden
	 * 			by the derived classes when I pass a potential_obj*
	 */
	virtual ~potential_obj();
	virtual double phi(double R,double z);
	virtual double dphidR(double R,double z);
	virtual double dphidz(double R,double z);
	virtual double PhiR(double R);
	virtual double dPhiR(double R);
	virtual double dPhi(double* x2,double *f2);
	double M_coeff();
	double b_coeff();
private:
	int whichpotential_;
	int npar_;
	char *fname_;
	GalaxyPotential *PhiWD_;
	double a_[5]; //Coefficients (M,b,q,etc)
};



#endif /* POTENTIAL_H_ */
