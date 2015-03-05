/*
 * potLeg.h
 *
 *  Created on: Feb 24, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_POTLEG_H_
#define INCLUDE_POTLEG_H_

#include "Potential.h"

void computeNewPhi(Potential *,double **rhlH=rhl,double ** sigRlH=sigRl,
		double **sigplH=sigpl,double ** sigzlH=sigzl,double **sigRzlH=sigRzl,
		double **vrotlH=vrotl);


#endif /* INCLUDE_POTLEG_H_ */
