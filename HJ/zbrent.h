/*
 * zbrent.h
 *
 *  Created on: Apr 3, 2012
 *      Author: fermani
 */

#ifndef ZBRENT_H_
#define ZBRENT_H_

#include "potential_obj.h"
#include "EtildeJ.h"

double zbrent_axR(double eq_sm,double at_this_R, double (*func)(double,double), double x1, double x2, double tol);

#endif /* ZBRENT_H_ */
