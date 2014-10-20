//for f(x,y,z) oct_int evaluates rho, or rho & rho*ybar, or rho, rho*ybar, rho*y^2bar or rho, rho*ybar, and all three rho*x_i^2bar 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
using namespace std;

#define MAX(A,B) ((A)>(B)?(A):(B))
#define SMALL 1.e-7

#include "uv_orb2.h"
#include "node1.h"
#include "node2.h"
#include "node3.h"

double oct_int(double (*fn)(double),double x0,double x1,int lmax){
	double I=0, DX=x1-x0;
	node1 domain(x0,0,DX);
	domain.set_values(fn);
	eval_node1(fn,domain,DX,&I,lmax);
	return I;
}
double oct_int(uv_orb *orb,double (uv_orb::*fn)(double),double x0,double x1,int lmax){
	double I=0, DX=x1-x0;
	node1 domain(x0,0,DX);
	domain.set_values(orb,fn);
	eval_node1(orb,fn,domain,DX,&I,lmax);
	return I;
}
double oct_int(double (*fn)(double,double),double x0,double x1,
	       double y0,double y1,int lmax){
	double I=0;
	location2 DIFF={x1-x0,y1-y0};
	node2 domain(x0,y0,0,DIFF);
	domain.set_values(fn);
	eval_node2(fn,domain,DIFF,&I,lmax);
	return I;
}
double oct_int(double (*fn)(double,double),double x0,double x1,
	       double y0,double y1,int lmax,double *J){
	double I=0; *J=0;
	location2 DIFF={x1-x0,y1-y0};
	node2 domain(x0,y0,0,DIFF);
	domain.set_values(fn);
	eval_node2(fn,domain,DIFF,&I,J,lmax);
	return I;
}
double oct_int(double*x,double (*fn)(double*,double,double),double x0,double x1,
	       double y0,double y1,int lmax){
	double I=0;
	location2 DIFF={x1-x0,y1-y0};
	node2 domain(x0,y0,0,DIFF);
	domain.set_values(x,fn);
	eval_node2(x,fn,domain,DIFF,&I,lmax);
	return I;
}
double oct_int(double*x,double (*fn)(double*,double,double),double x0,double x1,
	       double y0,double y1,int lmax,double *J){
	double I=0; *J=0;
	location2 DIFF={x1-x0,y1-y0};
	node2 domain(x0,y0,0,DIFF);
	domain.set_values(x,fn);
	eval_node2(x,fn,domain,DIFF,&I,J,lmax);
	return I;
}
double oct_int(double (*fn)(double,double,double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(fn);
	eval_node3(fn,domain,DIFF,&I,lmax);
	return I;
}
double oct_int(double *x,double (*fn)(double*,double,double,double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(x,fn);
	eval_node3(x,fn,domain,DIFF,&I,lmax);
	return I;
}
double oct_int(double (*fn)(double,double,double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax,double *J){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(fn);
	eval_node3(fn,domain,DIFF,&I,J,lmax);
	return I;
}
double oct_int(double *x,double (*fn)(double*,double,double,double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax,double *J){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(x,fn);
	eval_node3(x,fn,domain,DIFF,&I,J,lmax);
	return I;
}
double oct_int(double (*fn)(double,double,double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax,double *J,double *L){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0; *L=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(fn);
	eval_node3(fn,domain,DIFF,&I,J,L,lmax);
	return I;
}
double oct_int(double *x,double (*fn)(double*,double,double,double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax,double *J,double *L){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0; *L=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(x,fn);
	eval_node3(x,fn,domain,DIFF,&I,J,L,lmax);
	return I;
}
double oct_int(double (*fn)(double,double,double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax,double *J,double *K,double *L,double *M){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0; *K=0; *L=0; *M=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(fn);
	eval_node3(fn,domain,DIFF,&I,J,K,L,M,lmax);
	return I;
}
double oct_int(double *x,double (*fn)(double*,double,double,double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax,double *J,double *K,double *L,double *M){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0; *K=0; *L=0; *M=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(x,fn);
	eval_node3(x,fn,domain,DIFF,&I,J,K,L,M,lmax);
	return I;
}
double oct_int(double *x,double (*fn)(double*,double,double,double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax,
	       double *J,double *K,double *L,double *M,double *N){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0; *K=0; *L=0; *M=0; *N=0;
	node3 domain(x0,y0,z0,0,DIFF);
	domain.set_values(x,fn);
	eval_node3(x,fn,domain,DIFF,&I,J,K,L,M,N,lmax);
	return I;
}
double oct_int(double *x,double (*fn)(double*,double,double,double),
	       double (*fL)(double),double x0,double x1,
	       double y0,double y1,double z0,double z1,int lmax,
	       double *J,double *K,double *L,double *M,double *N){
	location3 DIFF={x1-x0,y1-y0,z1-z0};
	double I=0; *J=0; *K=0; *L=0; *M=0; *N=0;
	node3 domain(x0,y0,z0,0,DIFF);	domain.set_values(x,fn);
	eval_node3(x,fn,fL,domain,DIFF,&I,J,K,L,M,N,lmax);
	return I;
}
