/*
 * integ.h
 *
 *  Created on: Feb 23, 2015
 *      Author: L. Posti
 */

#ifndef INCLUDE_INTEG_H_
#define INCLUDE_INTEG_H_

#include "Potential.h"

double rhofDF(double, double, Potential *);
double rhofDF(double, double, Potential *, double *, double *, double *, double *, double *);


#endif /* INCLUDE_INTEG_H_ */
