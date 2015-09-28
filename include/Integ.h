/*
 * integ.h
 *
 *  Created on: Feb 23, 2015
 *      Author: L. Posti
 */

#ifndef INCLUDE_INTEG_H_
#define INCLUDE_INTEG_H_

#include "Potential.h"

double SigmaDF(double, double);
double rhofDF(double, double, Potential *);
double rhofDF(double, double, Potential *, double *, double *, double *, double *, double *);
double line_profile(double,double,double,Potential *);

#endif /* INCLUDE_INTEG_H_ */
