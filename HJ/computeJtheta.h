/*
 * computeJtheta.h
 *
 *  Created on: 27/giu/2014
 *      Author: morpheus
 */

#ifndef COMPUTEJTHETA_H_
#define COMPUTEJTHETA_H_

// private methods
//double computeJtheta(eveq_p * par, double Rc);

// public high-level functions
double getJtheta   (eveq_p * par, double Rc);
vec2   getJthetaJr (eveq_p * par, double Rc);


#endif /* COMPUTEJTHETA_H_ */
