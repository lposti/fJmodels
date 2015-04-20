/*
 * Uvpt.h
 *
 *  Created on: Feb 22, 2015
 *      Author: L. Posti
 */

#ifndef INCLUDE_UVPT_H_
#define INCLUDE_UVPT_H_

class pt{
public:
	/* Constructors */
	pt(const double, const double, const double);    // get (u,v,pu,pv) given (R,z,pR,pz)
	pt(const double, const double*, const double*);  // get (R,z) given (u,v)

	/* data */
	double Delta,Delta2;
	double R,z,u,v;
	double shu,chu,shu2,chu2;
	double sv,cv,sv2,cv2;
	double pR,pz,p2,pu,pv;

private:
	double R2,z2,r2;
};


#endif /* INCLUDE_UVPT_H_ */
