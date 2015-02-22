/*
 * Potential.h
 *
 *  Created on: Feb 15, 2015
 *      Author: morpheus
 */

#ifndef POT_POTENTIAL_H_
#define POT_POTENTIAL_H_

#include <iostream>
#include <cstring>

class Potential {

public:
	Potential();
	virtual ~Potential();

	/* data */
	double ** __restrict poly, ** __restrict I_int, ** __restrict I_ext;
	double * __restrict ci, * __restrict wi, * __restrict__ si, * __restrict pol;
	bool Ints,canEv;

	/* methods */
	// Guess density methods
	void selectGuessRho(const std::string);
	double (*rhoGuess)(double, double);
	void computeGuessRhl();
	// compute Phi and its derivative
	double operator () (double, double);
	void computePhil();
	double dr(const double,const double,double *);
	double dtheta(const double, const double);
	double dR(const double,const double);
	double dz(const double,const double);

private:
	void initLeg();
	void computeInts();

};

#endif /* POT_POTENTIAL_H_ */
