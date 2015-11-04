/*
 * GaussQuad.h
 *
 *  Created on: Feb 23, 2015
 *      Author: L. Posti
 */

#ifndef INCLUDE_GAUSSQUAD_H_
#define INCLUDE_GAUSSQUAD_H_

#include "GaussPts.h"

struct vec6d { double x0,x1,x2,x3,x4,x5; };

double Int1D_11     (double (*f)(double));
double Int1D_01     (double (*f)(double));
double Int1D_01     (double (*f)(double,void*), void *);
double Int2D_1111   (double (*f)(double,double));
double Int2D_0111   (double (*f)(double,double,void*), void *);
double Int2D_0101   (double (*f)(double,double,void*), void *);
double Int3D_011101 (double (*f)(double,double,double));
double Int3D_010101 (double (*f)(double,double,double,void*), void *);
double Int3D_011101 (double (*f)(double,double,double,void*), void *);
double Int3D_111101 (double (*f)(double,double,double,void*), void *);
double Int3D_111111 (double (*f)(double,double,double,void*), void *);
double Int3D_011101_vec (struct vec6d (*f)(double,double,double,void*), void *, double *);
double Int3D_111111_vec (struct vec6d (*f)(double,double,double,void*), void *, double *);



#endif /* INCLUDE_GAUSSQUAD_H_ */
