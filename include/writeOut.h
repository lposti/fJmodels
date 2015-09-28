/*
 * writeOut.h
 *
 *  Created on: Feb 27, 2015
 *      Author: L. Posti
 */

#ifndef INCLUDE_WRITEOUT_H_
#define INCLUDE_WRITEOUT_H_

#include <string>
#include <sstream>
#include <iomanip>
#include "Grid.h"
#include "Potential.h"
#include "readParam.h"

//#define LINEPROFILE

void writeOut(const struct fJParams&, const int,
		const unsigned comp=1, double **rhlH=rhl,double **sigRlH=sigRl,
		double **sigplH=sigpl,double **sigzlH=sigzl,double **sigRzlH=sigRzl,
		double **vrotlH=vrotl,double **philH=phil, double **PrH=Pr, double **Pr2H=Pr2);

void writeOut(const struct fJParams&, const int, Potential *p,
		const unsigned comp=1, double **rhlH=rhl,double **sigRlH=sigRl,
		double **sigplH=sigpl,double **sigzlH=sigzl,double **sigRzlH=sigRzl,
		double **vrotlH=vrotl,double **philH=phil, double **PrH=Pr, double **Pr2H=Pr2);


template <typename T> std::string toString( const T& t, int eps=2){

	std::stringstream ss;
	ss << std::fixed;				// format: could be also scientific
	ss << std::setprecision(eps);
	ss << t;
	return ss.str();
}


#endif /* INCLUDE_WRITEOUT_H_ */
