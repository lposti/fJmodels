/*
 * GaussQuad.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: morpheus
 */

#define QUADORD 12

static double xpt11[QUADORD]={-0.981560634,-0.904117256,-0.769902674,
                              -0.587317954,-0.367831499,-0.125233409,
                               0.125233409,0.367831499,0.587317954,
                               0.769902674,0.904117256,0.981560634},
              xw11[QUADORD]  ={0.047175336387,0.10693932600,0.16007832854,
                               0.20316742672,0.23349253654,0.24914704581,
                               0.24914704581,0.23349253654,0.20316742672,
                               0.16007832854,0.10693932600,0.047175336387};

static double xpt01[QUADORD] ={0.00921968288,0.0479413718,0.1150486629,
                               0.2063410229,0.316084251,0.437383296,
                               0.562616704,0.683915749,0.793658977,
                               0.884951337,0.952058628,0.990780317},
              xw01[QUADORD]  ={0.023587668193,0.053469662998,0.08003916427,
                               0.10158371336,0.11674626827,0.12457352291,
                               0.12457352291,0.11674626827,0.10158371336,
                               0.08003916427,0.053469662998,0.023587668193};

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
