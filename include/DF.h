/*
 * df.h
 *
 *  Created on: Feb 23, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_DF_H_
#define INCLUDE_DF_H_

#include "Potential.h"
#include "readParam.h"

void setDF(const double, const double, const double, const double, Potential *);
void setDF(struct fJParams, Potential *,const unsigned comp=1);
double df(const double *, const double *);




#endif /* INCLUDE_DF_H_ */
