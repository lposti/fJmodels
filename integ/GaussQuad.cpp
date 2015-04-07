/*
 * GaussQuad.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: morpheus
 */

#include "GaussQuad.h"
#include "GaussPts.h"

double Int1D_11 (double (*f)(double)){
	double sum=0.;

#pragma omp parallel for reduction(+:sum)
	for (unsigned i=0; i<QUADORD; i++)
		sum+= xw11[i]*f(xpt11[i]);

	return sum;
}

double Int2D_1111 (double (*f)(double,double)){
	double sum=0;

#pragma omp parallel for collapse(2) reduction(+:sum)
	for (unsigned i=0; i<QUADORD; i++)
		for (unsigned j=0; j<QUADORD; j++)
			sum+=xw11[i]*xw11[j]*f(xpt11[i],xpt11[j]);

	return sum;
}

double Int2D_0111 (double (*f)(double,double)){
	double sum=0;

#pragma omp parallel for collapse(2) reduction(+:sum)
	for (unsigned i=0; i<QUADORD; i++)
		for (unsigned j=0; j<QUADORD; j++)
			sum+=xw01[i]*xw11[j]*f(xpt01[i],xpt11[j]);
	return sum;
}

double Int3D_011101 (double (*f)(double,double,double)){
	double sum=0;

#pragma omp parallel for collapse(3) reduction(+:sum)
	for (unsigned i=0; i<QUADORD; i++)
		for (unsigned j=0; j<QUADORD; j++)
			for (unsigned k=0; k<QUADORD; k++)
				sum+=xw01[i]*xw11[j]*xw01[k]*f(xpt01[i],xpt11[j],xpt01[k]);

	return sum;
}

/* Version with paramters passed as void ptr. */
double Int3D_011101 (double (*f)(double,double,double,void*), void * params){
	double sum=0;

#pragma omp parallel for collapse(3) reduction(+:sum)
	for (unsigned i=0; i<QUADORD; i++)
		for (unsigned j=0; j<QUADORD; j++)
			for (unsigned k=0; k<QUADORD; k++)
				sum+=xw01[i]*xw11[j]*xw01[k]*f(xpt01[i],xpt11[j],xpt01[k],params);

	return sum;
}

/* Version with also vrot and sigma integration: a double ptr with dimension 6 must be passed */
double Int3D_011101_vec (struct vec6d (*f)(double,double,double,void*), void * params, double * out){
	double sum=0,sum1=0,sum2=0,sum3=0,sum4=0,sum5=0;
	struct vec6d res;

#pragma omp parallel for collapse(3) reduction(+:sum,sum1,sum2,sum3,sum4,sum5) private (res)
	for (unsigned i=0; i<QUADORD; i++)
		for (unsigned j=0; j<QUADORD; j++)
			for (unsigned k=0; k<QUADORD; k++){
				res=f(xpt01[i],xpt11[j],xpt01[k],params);
				sum +=xw01[i]*xw11[j]*xw01[k]*res.x0;
				sum1+=xw01[i]*xw11[j]*xw01[k]*res.x1;
				sum2+=xw01[i]*xw11[j]*xw01[k]*res.x2;
				sum3+=xw01[i]*xw11[j]*xw01[k]*res.x3;
				sum4+=xw01[i]*xw11[j]*xw01[k]*res.x4;
				sum5+=xw01[i]*xw11[j]*xw01[k]*res.x5;
			}

	out[0]=sum;  out[1]=sum1; out[2]=sum2;
	out[3]=sum3; out[4]=sum4; out[5]=sum5;
	return sum;
}
