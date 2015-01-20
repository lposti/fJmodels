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
//#define JAFFE

#define EXTERNPOT
#define CONSTOMRATIO
#define HJGJ

#define TPI 6.283185307179586
extern double *EofJ,*LzofJ;
extern int rep;

extern double mass,J0,r0;

void setdf(double dr_in, double dphi_in, double dz_in);
void setdf(double dr_in, double dphi_in, double dz_in,
		double dr_g_in, double dphi_g_in, double dz_g_in);
void setMJ0(double M_in, double J0_in);

double df(double *x,double *v);


#endif /* ISODF4_H_ */
