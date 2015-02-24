/*
 * GaussQuad.h
 *
 *  Created on: Feb 23, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_GAUSSQUAD_H_
#define INCLUDE_GAUSSQUAD_H_

double Int1D_11     (double (*f)(double));
double Int2D_1111   (double (*f)(double,double));
double Int2D_0111   (double (*f)(double,double));
double Int3D_011101 (double (*f)(double,double,double));
double Int3D_011101 (double (*f)(double,double,double,void*), void *);



#endif /* INCLUDE_GAUSSQUAD_H_ */
