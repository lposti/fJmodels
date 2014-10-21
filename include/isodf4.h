/*
 * isodf4.h
 *
 *  Created on: 05/feb/2014
 *      Author: morpheus
 */

#ifndef ISODF4_H_
#define ISODF4_H_

//#define ISOCHRONE
//#define NFW
#define HERNQUIST
//#define ISOTHERMAL

#define CONSTOMRATIO
#define HJGJ

extern double *EofJ,*LzofJ;
extern int rep;

void setdf(double dr_in, double dphi_in, double dz_in);

double df(double *x,double *v);


#endif /* ISODF4_H_ */
