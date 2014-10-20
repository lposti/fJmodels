/*
 * int_eq.h
 *
 *  Created on: Mar 12, 2013
 *      Author: Fermani
 */

#ifndef INT_EQ_H_
#define INT_EQ_H_

void integrate_eq(int T,double **orbit,int dim,int (*eq)(double,const double*, double*,void*),double *pt,void *par, double eps);

void integrate_eq_lastpt(int T,double *lastpt,int dim,int (*eq)(double,const double*, double*,void*),double *pt,void *par, double eps);

double trace(int dim, double **matrix);

double determinant33(double **matrix);

void multiplic_scalar(int dim, double scalar, double **matrix); // sores new values in matrix

void sum_matrices(int dim,double **matrix1,double **matrix2);

void subtract_matrices(int dim,double **matrix1,double **matrix2);

void eigenvalues_33(double *eig1, double *eig2, double *eig3, double **M);

#endif /* INT_EQ_H_ */
