/*
 * tools.cpp
 *
 *  Created on: Mar 12, 2013
 *      Author: Fermani
 */

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "tools.h"
#include "Pi.h"

using namespace std;

using std::cout;

void integrate_eq(int T, double **orbit, int dim,
		int (*eq)(double, const double*, double*, void*), double *pt, void *par,
		double eps) {
	size_t N = dim;
	gsl_odeiv2_system sys = { eq, NULL, N, par };
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,
			gsl_odeiv2_step_rk8pd, eps, eps, 0.0);

	int i, j;
	double t;
	double t0 = 0.0;
	double y[dim];
	for (i = 0; i < dim; i++)
		y[i] = pt[i];

	int status;
	for (i = 0; i != T; i++) {
		t = double(i);
		status = gsl_odeiv2_driver_apply(d, &t0, t, y);
		if (status != GSL_SUCCESS) {
			cout << "error, return value= " << status << "\n";
			exit(1);
		}
		for (j = 0; j < dim; j++)
			orbit[i][j] = y[j];
	}

	gsl_odeiv2_driver_free(d);
}

void integrate_eq_lastpt(int T, double *lastpt, int dim,
		int (*eq)(double, const double*, double*, void*), double *pt, void *par,
		double eps) {
	size_t N = dim;
	gsl_odeiv2_system sys = { eq, NULL, N, par };
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,
			gsl_odeiv2_step_rk8pd, eps, eps, 0.0);

	int i, j;
	double t;
	double t0 = 0.0;
	double y[dim];
	for (i = 0; i < dim; i++)
		y[i] = pt[i];

	int status;
	for (i = 0; i != T; i++) {
		t = double(i);
		status = gsl_odeiv2_driver_apply(d, &t0, t, y);
		if (status != GSL_SUCCESS) {
			cout << "error, return value= " << status << "\n";
			exit(1);
		}
		for (j = 0; j < dim; j++)
			lastpt[j] = y[j];
	}

	gsl_odeiv2_driver_free(d);
}

double trace(int dim,double **matrix){
	double res=0;
	for(int i=0;i<dim;i++) res+=matrix[i][i];
	return res;
}

double determinant33(double **matrix){
	double p1,p2,p3,m1,m2,m3;
	p1=matrix[0][0]*matrix[1][1]*matrix[2][2];
	p2=matrix[0][1]*matrix[1][2]*matrix[2][0];
	p3=matrix[0][2]*matrix[1][0]*matrix[2][1];

	m1=matrix[2][0]*matrix[1][1]*matrix[0][2];
	m2=matrix[2][1]*matrix[1][2]*matrix[0][0];
	m3=matrix[2][2]*matrix[1][0]*matrix[0][1];
	return p1+p2+p3-(m1+m2+m3);
}

void multiplic_scalar(int dim,double scalar,double **matrix){
	int i,j;
	for(i=0;i<dim;i++) for(j=0;j<dim;j++) matrix[i][j]*=scalar;
}

void sum_matrices(int dim,double **matrix1,double **matrix2){
	int i,j;
	for(i=0;i<dim;i++) for(j=0;j<dim;j++) matrix1[i][j]=matrix1[i][j]+matrix2[i][j];
}

void subtract_matrices(int dim,double **matrix1,double **matrix2){
	int i,j;
	for(i=0;i<dim;i++) for(j=0;j<dim;j++) matrix1[i][j]=matrix1[i][j]-matrix2[i][j];
}

void eigenvalues_33(double *eig1, double *eig2, double *eig3, double **M){
	double p,q,r,phi; // from wiki-code
	int i,j;
	double **I= new double* [3];
	for(i=0;i<3;i++) I[i]=new double[3];
	for(i=0;i<3;i++) for(j=0;j<3;j++) I[i][j]=0.; //Initialize identify matrix
	for(j=0;j<3;j++) I[j][j]=1.0;

	p = pow(M[0][1],2.0) + pow(M[0][2],2.0) + pow(M[1][2],2.0);
	//cout <<  "p " << p << "\n";
	if(p == 0){
		// A is diagonal.
		*eig1 = M[0][0];
		*eig2 = M[1][1];
		*eig3 = M[2][2];
	}
	else{
		q=trace(3,M)/3.;
		p=pow(M[0][0]-q,2.0)+pow(M[1][1]-q,2.0)+pow(M[2][2]-q,2.0)+2.0*p;
		p=sqrt(p/6.);
		multiplic_scalar(3,-q,I);
		sum_matrices(3,M,I);
		multiplic_scalar(3,1.0/p,M);
		r=determinant33(M)/2.0;

		//cout <<  "r " << r << "\n";

		// In exact arithmetic for a symmetric matrix  -1 <= r <= 1
		// but computation error can leave it slightly outside this range.
		if (r<=-1) phi = Pi/3.;
		else{
			if(r>=1) phi = 0.;
			else phi = acos(r)/3.;
		}
		//cout <<  "phi " << phi << "\n";
		//the eigenvalues satisfy eig3 <= eig2 <= eig1
		*eig1 = q + 2. * p * cos(phi);
		*eig3 = q + 2. * p * cos(phi + Pi * (2./3.));
		*eig2 = 3. * q - *eig1 - *eig3;     // since trace(A) = eig1 + eig2 + eig3
		//cout << "Eigenvalues: " << *eig1 << " " << *eig2 << " " << *eig3 << "\n";
	}
}
