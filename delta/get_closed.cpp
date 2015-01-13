#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "potLEG.h"
#include "press.h"
#include "mongo.h"

double get_delta2(double,double*,double*,int);

double dPhi_eff(double *par,double x){
	double x2[2]={x,0},dP[2];
	dPhi(x2,dP);
	return dP[0]-par[1]/pow(x,3);
}
void derivs(double *par,double t,double *y,double *dydt){
	double dP[2];
	dPhi(y,dP);
	dP[0]-=par[1]/pow(y[0],3);
	for(int i=0;i<2;i++){
		dydt[i]=y[i+2];
		dydt[i+2]=-dP[i];
	}
}
double vsqx(double *par,double x){
	return 2*(par[0]-(PhiR(x)+.5*par[1]/(x*x)));
}
double Ez(double *par,double *x2){
	return .5*(pow(x2[2],2)+pow(x2[3],2))+.5*par[1]/pow(x2[0],2)+Phi(x2[0],x2[1]);
}
double goround(double *par,double x){//throws up and determines Rdown
	double vsq=vsqx(par,x);
	if(vsq<0){
		printf("%f %f\n",x,vsq); exit(0);
	}
	double x2[4]={x,0,0,sqrt(vsq)};
	double dxdt[4],xscal[4],t,hdid,hnext,htry=2.e-2,eps=1.e-10;
	t=0;
	for(int i=0; i<2; i++){
		xscal[i]=x; xscal[i+2]=fabs(x2[3]);
	}
	double Rlast,zlast;
	while(x2[1]>=0){
		Rlast=x2[0]; zlast=x2[1];
		derivs(par,t,x2,dxdt);
		rkqs(par,x2,dxdt,4,&t,htry,eps,xscal,&hdid,&hnext,&derivs); htry=hnext; 
	}
	return x-(x2[0]-x2[1]/(zlast-x2[1])*(Rlast-x2[0]));//difference between up and down
}
int go_up(double *par,double x,double *Ri,double *zi,int nmaxR){//throws up and determines Rdown
	double vsq=vsqx(par,x);
	if(vsq<0){
		printf("%f %f\n",x,vsq); exit(0);
	}
	double x2[4]={x,0,0,sqrt(vsq)};
	double dxdt[4],xscal[4],t,hdid,hnext,htry=2.e-2,eps=1.e-10;
	t=0;
	for(int i=0; i<2; i++){
		xscal[i]=x; xscal[i+2]=fabs(x2[3]);
	}
	int j=0;
	while(x2[3]>=0){
		derivs(par,t,x2,dxdt);
		rkqs(par,x2,dxdt,4,&t,htry,eps,xscal,&hdid,&hnext,&derivs); htry=hnext;
		Ri[j]=x2[0]; zi[j]=x2[1]; j++; if(j==nmaxR) break;
	}
	//if(Ri[0]<zi[j-1]) setcolour("black"); else setcolour("red");
	return j;
}
void plot_updown(double *par,double x){//throws up and determines Rdown
	double vsq=vsqx(par,x),Ri[1000],zi[1000];
	double x2[4]={x,0,0,sqrt(vsq)};
	double dxdt[4],xscal[4],t,hdid,hnext,htry=2.e-2,eps=1.e-10;
	t=0;
	for(int i=0; i<2; i++){
		xscal[i]=x; xscal[i+2]=fabs(x2[3]);
	}
	int j=0;
	while(x2[1]>=0){
		Ri[j]=x2[0]; zi[j]=x2[1]; j++;
		derivs(par,t,x2,dxdt);
		rkqs(par,x2,dxdt,4,&t,htry,eps,xscal,&hdid,&hnext,&derivs); htry=hnext; 
	}
	connect(Ri,zi,j);
}
double sorted(double *par,double x0){
	int k=0;
	while(vsqx(par,x0)<0){
		if(dPhi_eff(par,x0)<0) x0*=1.1;
		else x0*=.9;
		k++; if(k>50){ printf("vsq: %f %f %f\n",x0,par[0],par[1]); exit(0);}
	}
	return x0;
}
#define nmaxRi 200
double get_closed(double x0,double Lz,double E1){
	double par[2]={E1,Lz*Lz};
	x0=sorted(par,x0);// ensure vsq>=0
	double x1,dx0,dx1;
	dx0=goround(par,x0);
	double dx=.8; x1=sorted(par,x0/dx); dx1=goround(par,x1);
	int k=0,bo=0;
	while(dx0*dx1>0){//we haven't bracketed the root
		if(fabs(dx1)<fabs(dx0)){//more of same
			x0=x1; dx0=dx1;
		}else{// back off
			if(bo==0){
				dx=1/dx; bo=1;
			}else dx=pow(dx,-0.9);
		}
		x1=sorted(par,x0*dx); dx1=goround(par,x1);
		k++;
		if(k==21){
			printf("get_closed: %f %f %g %g %g\n",x0,x1,dx0,dx1,dx);
			if(fabs(dx0)<1.e-7){
				return -x0;
			}else{
				printf("problem in get_closed(): %f %f %f %f\n",x0,x1,dx0,dx1);
				return -1;
			}
		}
	}
	double R0=zbrent(par,&goround,x0,x1,dx0,dx1,1.e-6,25);
	double Ri[nmaxRi],zi[nmaxRi];
	int np=go_up(par,R0,Ri,zi,nmaxRi);
	//plot_updown(par,R0);
	return get_delta2(R0,Ri,zi,np);
}
