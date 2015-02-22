/*
 * Delta.h
 *
 *  Created on: Feb 16, 2015
 *      Author: morpheus
 */

#ifndef DELTA_H_
#define DELTA_H_

#include "Potential.h"

/* Parameters to be passed to Delta */
struct eqsPar { Potential p; double E,L,Lz,R0;};

/*
 *  Interface class for passing "condition functor" to
 *  evolveOrbit method.
 */
class condInterf {
public:
	virtual ~condInterf() {}
	virtual bool operator() (double *,int,int) = 0;
};

/* Non-abstract class for condition z>=0 */
class Cond_zgt0 : condInterf {
public:
	virtual bool operator() (double y[], int i, int size) {
		if (i<size-1 and y[1]>=0.) return true;
		else return false;
	}
};

/* Non-abstract class for condition pz>=0 */
class Cond_pzgt0 : condInterf {
public:
	virtual bool operator() (double y[], int i, int size) {
		if (i<size-1 and y[3]>=0.) return true;
		else return false;
	}
};

/*
 *  Delta class: computes best focal distance for ellipsoidal coordinate
 *  system for a given Energy and Angular momentum.
 */
class Delta {
public:

	// Constructor & Destructor
	Delta(struct eqsPar * par_in);
	virtual ~Delta() {}

	// data
	struct eqsPar * par;
	double R0;

	// methods
	inline double RcEfn(double);
	//double RcE(double);
	inline double ds2dv (const double);
	double ds2dD2TOT(double, void *);
	unsigned evolveOrbit(double,double *,double *,double *,double *,
			             const unsigned,condInterf *);
	double getDR(double);
	double getBestDelta();

private:

	// data
	double Rc,D2,sqDR,v0,cv,sv,Ri,zi;
	unsigned size;

	// methods
	double getv();
	double getJr0();
	double ensureVsq(double);
	double dPhiEff(const double);
};

/* Parameters needed by GSL functions */
struct Rzptr { double *R,*z; int N; Delta * dd;};


#endif /* DELTA_H_ */
