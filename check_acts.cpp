#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "potLEG.h"
#include "press.h"
#include "uv_orb3.h"
#include "tables3.h"

void derivs(double*,double,double*,double*);
void check_acts(double *x,double *v){//integrates from given point & determins rms of actions
	double x2[4]={x[0],x[1],v[0],v[2]},Lz=x[0]*v[1];
	double dxdt[4],xscal[4],hdid,hnext,htry=2.e-2,eps=1.e-10;
	double H=.5*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2))+Phi(x[0],x[1]);
	double par[2]={H,Lz*Lz};
	double r=sqrt(pow(x[0],2)+pow(x[1],2));
	double t=0,Jrbar=0,Jzbar=0,Jrsig=0,Jzsig=0,Jre=0,Jze=0;
	for(int i=0; i<2; i++){
		xscal[i]=r; xscal[i+2]=r;
	}
	int j=0,step=10,nt=201;
	while(j<nt){
		derivs(par,t,x2,dxdt);
		rkqs(par,x2,dxdt,4,&t,htry,eps,xscal,&hdid,&hnext,&derivs);
		htry=hnext;
		double p[2]={x2[2],x2[3]},Jr,Jz,Jr1,Jz1;
		uv_orb uvorb(Deltafn(H),Lz,Phi(x2[0],x2[1]),x2,p);
		int nt=int_acts(Lz,uvorb.E,uvorb.Er,&Jr,&Jz);
		Jr1=uvorb.Ju(); Jz1=uvorb.Jv();
		Jrbar+=Jr1; Jzbar+=Jz1;	Jrsig+=pow(Jr1,2); Jzsig+=pow(Jz1,2);
		if(nt==0){
			if(j==0) printf("off grid ");
		}else{
			Jre+=pow(Jr1-Jr,2); Jze+=pow(Jz1-Jz,2);
			if(j==0) printf("(%f %f) (%f %f)\n",Jr,Jr1,Jz,Jz1);
		}
		j++;
	}
	Jrbar/=nt; Jzbar/=nt;
	Jrsig=sqrt(Jrsig/nt-Jrbar*Jrbar)/(Jrbar+Jzbar);
	Jzsig=sqrt(Jzsig/nt-Jzbar*Jzbar)/(Jrbar+Jzbar);
	Jre/=nt; Jre=sqrt(Jre)/(Jrbar+Jzbar); Jze/=nt; Jze=sqrt(Jze)/(Jrbar+Jzbar);
	printf("From (%f %f) to (%f %f)\n",x[0],x[1],x2[0],x2[1]);
	printf("Jb: (%f %f) Js: (%g %g) Je: (%f %f)\n",Jrbar,Jzbar,Jrsig,Jzsig,Jre,Jze); 
}
