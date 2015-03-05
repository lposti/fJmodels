/*
 * models.h
 *
 *  Created on: Feb 15, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_MODELS_H_
#define INCLUDE_MODELS_H_

#include "readParam.h"

// q=flattening, r0=break radius, J0=scale action
extern double q,q2;
extern double mass,r0,J0;

void SetModel(struct fJParams, const unsigned comp=1);
double rhoHern(double, double);
double rhoIsoch(double, double);
double rhoNFW(double, double);
double rhoNFWext(double, double);



#endif /* INCLUDE_MODELS_H_ */
