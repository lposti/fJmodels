#include <stdio.h>
#include <math.h>
#include "potLEG.h"
#include "mongo.h"
#include "press.h"

static double *Rp,*zp,DS2last;
static int N;
static double R,z,v0,sv,cv,R0,D2,sqDR;
static const double PIH=.5*acos(-1);
static double s2(double v){
	return pow(R0*sin(v)-R,2)+pow(sqDR*cos(v)-z,2);
}
static double ds2dv0(double v){
	return 2*(z*sqDR*sin(v)-cos(v)*(R*R0+D2*sin(v)));
}
static double get_v0(double Ri,double zi){//returns v of closest approach to (R,z)
	R=Ri; z=zi;
	if(ds2dv0(0)*ds2dv0(PIH)>0){
		printf("In get_v0 %f %f %f %f %f",Ri,zi,sqDR,R,R0);
	}
	v0=zbrent(&ds2dv0,0,PIH,ds2dv0(0),ds2dv0(PIH),1.e-10,25);
	sv=sin(v0); cv=cos(v0);
	return v0;
}
static double dv0dD2(void){
	return .5*sv*(z/sqDR-2*cv)/(cv*z*sqDR+sv*R*R0-D2*(cv*cv-sv*sv));
}
static double ds2dD2(void){
	return (cv-z/sqDR)*cv+ds2dv0(v0)*dv0dD2();
}
/*static double s2total(double D20){
	D2=D20; sqDR=sqrt(D2+R0*R0);
	double sum=0;
	for(int i=0; i<N; i++){
		get_v0(Rp[i],zp[i]); 
		sum+=s2(v0);
	}
	return sum;
}*/
static double ds2totaldD2(double D20){
	D2=D20; sqDR=sqrt(D2+R0*R0);
	/*testing stuff
	double v=.3,eps=.001; printf("%f ",ds2dv0(v)); 
	double ds=s2(v+eps); ds-=s2(v-eps); ds*=.5/eps; printf("%f\n",ds);
	get_v0(Rp[0],zp[0]); printf("%f ",ds2dD2());
	D2=D20+eps; sqDR=sqrt(D2+R0*R0); get_v0(Rp[0],zp[0]); ds=s2(v0);
	D2=D20-eps; sqDR=sqrt(D2+R0*R0); get_v0(Rp[0],zp[0]); ds-=s2(v0);
	ds*=.5/eps; printf("%f\n",ds);
	//end tests*/
	double sum=0;
	for(int i=0; i<N; i++){
		get_v0(Rp[i],zp[i]); 
		sum+=ds2dD2();
	}
	return sum;
}
#define NP 50
double get_delta2(double R,double *Ri,double *zi,int ni){
	R0=R; Rp = new double[ni]; zp = new double[ni]; N=ni;
	for(int i=0;i<ni; i++){
		Rp[i]=Ri[i]; zp[i]=zi[i];
	}
	double D2min=-.001,D2max=5,dsdDmin=ds2totaldD2(D2min),dsdDmax=ds2totaldD2(D2max);
/*	if(Rp[0]>zp[ni-1]){
		printf("(%f %f %g %g) ",Rp[0],zp[ni-1],Phi(Rp[0],zp[0]),Phi(Rp[ni-1],zp[ni-1]));
		float x[NP],y[NP],z[NP];
		D2=D2min;D2=DS2last;
		sqDR=sqrt(D2+R0*R0);
		for(int i=0; i<NP; i++){
			double v=i*PIH/(double)(NP-1);
			x[i]=R0*sin(v);	y[i]=sqDR*cos(v);
		}
		//plots(1,x,y,0,30,0,30,"R",1,"z",1,-.9,10);
		setcolour("green");
		connect(Ri,zi,ni);
		//grend(1);
	}*/
	while(dsdDmin*dsdDmax>0){
		if(dsdDmax<0){
			D2max*=5; dsdDmax=ds2totaldD2(D2max);
		}else{
			D2min-=.1; dsdDmin=ds2totaldD2(D2min);
		}
	}
	if(dsdDmin*dsdDmax<0)
		DS2last=D2=zbrent(&ds2totaldD2,D2min,D2max,dsdDmin,dsdDmax,1.e-10,25);

	if(Rp[0]>zp[ni-1] && DS2last<0.){
			FILE* forb=fopen("orb.dat","w");
			printf("(%f %f %g %g) ",Rp[0],zp[ni-1],Phi(Rp[0],zp[0]),Phi(Rp[ni-1],zp[ni-1]));
			float x[NP],y[NP],z[NP];
			D2=D2min;D2=DS2last;
			sqDR=sqrt(D2+R0*R0);
			for(int i=0; i<NP; i++){
				double v=i*PIH/(double)(NP-1);
				x[i]=R0*sin(v);	y[i]=sqDR*cos(v);
				fprintf(forb,"%f %f\n",x[i],y[i]);
			}
			//plots(1,x,y,0,30,0,30,"R",1,"z",1,-.9,10);
			//setcolour("green");
			//connect(Ri,zi,ni);
			printf(">> Delta=%f\n",DS2last);
			fclose(forb);

		}
	delete[] Rp; delete[] zp;
	return (DS2last);
}

/*
#define NP 50
int main(void){
	float x[NP],y[NP],z[NP];
	R0=.3; double D2min=.1,D2max=2;
	for(int i=0; i<NP; i++){
		x[i]=D2=D2min+i*D2max/(float)(NP-1);
		y[i]=ds2totaldD2(D2); z[i]=s2total(D2);
		//printf("%f %f %f\n",x[i],y[i],z[i]);
	}
	plots(NP,x,y,0,D2max,-0.4,0.4,"D",1,"ds",2,-.9,10);
	setcolour("red"); connect(x,z,NP); grend(1);
	D2=D2min;
	D2=zbrent(&ds2totaldD2,D2min,D2max,ds2totaldD2(D2min),ds2totaldD2(D2max),1.e-5,25);
	sqDR=sqrt(D2+R0*R0);
	printf("%f %f\n",D2,s2total(D2));
	for(int i=0; i<NP; i++){
		float v=i*PIH/(float)(NP-1);
		x[i]=R0*sin(v); y[i]=sqDR*cos(v);
	}
	plots(NP,x,y,0,1,0,1,"R",1,"z",1,-.9,10);
	points(43.15,1,Rp,zp,N); relocate(R0,0.01); point(33.15);
	grend(1);
}*/
