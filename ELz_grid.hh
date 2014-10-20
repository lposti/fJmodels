/*
 * ELz_grid.cpp
 *
 *  Created on: 04/feb/2014
 *      Author: LP
 */
#include "grid3.h"
#include "potLEG.h"
#include "ini_potLEG.h"

#define TINY 1.e-10
#define SMALLEST 10
#define MAX(A,B) ((A)>(B)?(A):(B))
#define SIGN(A,B) ((B)>=0 ?fabs(A):-fabs(A))
#define EPS 3.0e-8

double get_closed(double,double,double);

class ELz_grid
{
public:
	ELz_grid(double Rmax=50.);					// Constructor
	~ELz_grid();								// Destructor
	ELz_grid(const ELz_grid &L);             	// copy constructor
	ELz_grid & operator=(const ELz_grid &L); 	// assignment

	void make_NGRIDgrids();						// Makes Egrid, Dgrid
	void make_NRgrids();						// Makes Lzgrid, Egrid, Ecgrid, Vmax

	double *Lzgrid;								// Pointer to array Lz-grid
	double *Egrid;								// Pointer to array E-grid
	double *Ecgrid;								// Pointer to array E circ-grid
	double *Dgrid;								// Pointer to array Delta-grid
	double *Rgrid;								// Pointer to array Radial-grid
	double *Vmax;
	int gE,gR;									// Size of the grids
	double GetRc(double Lzsq,double R);
	double GetRc_ini(double Lzsq,double R);		// Rc for initial potential

private:
	int gEfn(int nL);
	double Rmin,R_max,Lzmin,Lzmax;
	double Rcfn(double *Lp,double R);
	double Rcfn_ini(double *Lp,double R);
	double RcEfn(double *Ep,double R);
	double RcE(double E,double Rc);
	double zbrent(double *pars,
	double (ELz_grid::*func)(double*,double),
	double x1, double x2,double fa,double fb,
	double tol,int itmax);
};


ELz_grid::ELz_grid(double Rmax){

	R_max=Rmax; Rmin=1e-4*R_max;
	Lzmax=R_max*0.98*sqrt(-2*PhiR(R_max));
	Lzmin=Rmin*sqrt(Rmin*dPhiR(Rmin));
	gR=NR; 	gE=NGRID;

	/* Array allocations */
	this->Rgrid = dmatrix(NR); this->Lzgrid = dmatrix(NR);
	this->Vmax = dmatrix(NR); this->Ecgrid = dmatrix(NR);
	this->Egrid = dmatrix(NGRID); this->Dgrid = dmatrix(NGRID);

}

void ELz_grid::make_NRgrids(){
	double R=.01*R_max;

	for(int i=0; i<gR; i++){
			this->Lzgrid[i]=Lzmin+i*(Lzmax-Lzmin)/(double)(gR-1);
			double Lzsq=pow(this->Lzgrid[i],2); R=this->Rgrid[i]=GetRc(Lzsq,R);
			double vc2=R*dPhiR(R);
			this->Ecgrid[i]=.5*vc2+PhiR(R);//Ecirc at R
			this->Vmax[i]=.98*sqrt(-2*this->Ecgrid[i]);//Range of Egrid at R
	}
}

void ELz_grid::make_NGRIDgrids(){

	double Lzgrid_h=Lzmin+(Lzmax-Lzmin)/(double)(gR-1);
	double Lzsq=pow(Lzgrid_h,2), R=GetRc(Lzsq,.01*R_max);
	double Ecgrid_h=.5*R*dPhiR(R)+PhiR(R);
	double Vmax_h=.98*sqrt(-2*Ecgrid_h);

	for(int nE=0; nE<gE; nE++){//loop over E
		double E=Ecgrid_h+.5*pow(Vmax_h*(nE+1)/(double)gE,2);
		this->Dgrid[nE]=get_closed(RcE(E,R),.25*Lzgrid_h,E);
		this->Egrid[nE]=E;
	}
}

ELz_grid::~ELz_grid(){
	delete [] Lzgrid; delete [] Ecgrid; delete [] Egrid;
	delete [] Dgrid; delete [] Rgrid; delete [] Vmax;
}

double ELz_grid::Rcfn(double *Lp,double R){//what vanishes at Rc for given Lz
	double x[2]={R,0},f[2]; dPhi(x,f);
	return (*Lp)/pow(R,3)-f[0];
}

double ELz_grid::Rcfn_ini(double *Lp,double R){//what vanishes at Rc for given Lz
	double x[2]={R,0},f[2]; dPhi_ini(x,f);
	return (*Lp)/pow(R,3)-f[0];
}

double ELz_grid::zbrent(double *pars,double (ELz_grid::*func)(double*,double),
						double x1, double x2,double fa,double fb,
						double tol,int itmax){
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
		printf("Root must be bracketed in zbrent:\nx1,x2,f1,f2: %g %g %g %g",x1,x2,fa,fb);
		exit(0);}
	fc=fb;
	for (iter=1;iter<=itmax;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a; fc=fa; e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b; b=c; c=a; fa=fb; fb=fc; fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s; q=1.0-s;
			} else {
				q=fa/fc; r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d; d=p/q;
			} else {
				d=xm; e=d;
			}
		} else {
			d=xm; e=d;
		}
		a=b; fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(this->*func)(pars,b);
	}
	if(fabs(xm)>100*tol1){
		printf("Too many iters in zbrent: %g %g ",b,fb);// exit(0);
	}
	return b;
}

double ELz_grid::GetRc(double Lzsq,double R){
	double Ri,Ro,dfRi,dfRo,dfRt=Rcfn(&Lzsq,R);
	if(dfRt>0){//inside
		while(dfRt>0){
			Ri=R; dfRi=dfRt; R*=1.2; dfRt=Rcfn(&Lzsq,R);
		} Ro=R; dfRo=dfRt;
	}else{//outside
		while(dfRt<0){
			Ro=R; dfRo=dfRt; R*=.8; dfRt=Rcfn(&Lzsq,R);
		} Ri=R; dfRi=dfRt;
	}
	return zbrent(&Lzsq,&ELz_grid::Rcfn,Ri,Ro,dfRi,dfRo,TINY,20);
}

double ELz_grid::GetRc_ini(double Lzsq,double R){
	double Ri,Ro,dfRi,dfRo,dfRt=Rcfn_ini(&Lzsq,R);
	if(dfRt>0){//inside
		while(dfRt>0){
			Ri=R; dfRi=dfRt; R*=1.2; dfRt=Rcfn_ini(&Lzsq,R);
		} Ro=R; dfRo=dfRt;
	}else{//outside
		while(dfRt<0){
			Ro=R; dfRo=dfRt; R*=.8; dfRt=Rcfn_ini(&Lzsq,R);
		} Ri=R; dfRi=dfRt;
	}
	return zbrent(&Lzsq,&ELz_grid::Rcfn_ini,Ri,Ro,dfRi,dfRo,TINY,20);
}

double ELz_grid::RcEfn(double *Ep,double R){//what vanishes at Rc for given E
	return .5*R*dPhiR(R)+PhiR(R)-(*Ep);
}

double ELz_grid::RcE(double E,double Rc){
	double R0=.8*Rc,f0=RcEfn(&E,R0);
	while(f0>0){
		R0*=.8; f0=RcEfn(&E,R0);
	}
	double R1=1.2*Rc, f1=RcEfn(&E,R1);
	while(f1<0){
		R1*=1.2;  f1=RcEfn(&E,R1);
	}
	return zbrent(&E,&ELz_grid::RcEfn,R0,R1,f0,f1,1.e-5,25);
}

int ELz_grid::gEfn(int nL){
	return MAX(SMALLEST,NGRID*(float)(NR-nL)/NR);//more E values at low Lz
}

