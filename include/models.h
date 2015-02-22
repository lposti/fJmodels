/*
 * models.h
 *
 *  Created on: Feb 15, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_MODELS_H_
#define INCLUDE_MODELS_H_

// q=flattening, r0=break radius, J0=scale action
extern const double q,q2;
extern const double mass,r0,J0;


double rhoHern(double, double);
double rhoIsoch(double, double);



#endif /* INCLUDE_MODELS_H_ */
