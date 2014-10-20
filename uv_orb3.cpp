#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "press.h"
#include "Rzuv.h"
#include "uv_orb3.h"
#include "potLEG.h"
#include "oct_int_exp.h"

#define TINY 1.e-12
#define SMALL 1.e-4
#define TPI 6.2831853072
#define PI  3.141592654
#define PIH 1.570796327
#define MAX(A,B) ((A)>(B)?(A):(B))
#define SIGN(A,B) ((B)>=0 ?fabs(A):-fabs(A))

static const double us=58,rtwo=sqrt(2); 
double oct_int(uv_orb*,double (uv_orb::*)(double),double,double,int);

#define EPS 3.0e-8

double Xbrent(uv_orb *orb,double (uv_orb::*func)(double), double x1, double x2,double fa,double fb,double tol,int itmax){
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
		printf("Root must be bracketed in Xbrent:\nx1,x2,f1,f2: %g %g %g %g",x1,x2,fa,fb);
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
		fb=(orb->*func)(b);
	}
	if(fabs(xm)>100*tol1){
		//printf("Too many iters in Xbrent: %g %g ",b,fb);// exit(0);
	}
	return b;
}

#undef EPS
double uv_orb::Phiu(double u){
	pt X(Delta,u,v);
	return Phi(X.R,X.z);
}
double uv_orb::dU(double u){//returns U(u)-U(u0)
	return (pow(sinh(u),2)+sv2)*Phiu(u)-(sh1sq+sv2)*Phiu1;
}
double uv_orb::dPhiu(double u,double *f){//returns Phiu(u) and dPhi/dR etc
	pt X(Delta,u,v);
	double x[2]={X.R,X.z};
	return dPhi(x,f);
}
double uv_orb::uturnfn(double u){//what vanishes at peri / apo
	double shu2=pow(sinh(u),2);
//	Iu=E*(X.shu2-sh1sq)-.5*pow(X.pu,2)/(X.Delta2)-.5*Lzsq/(X.Delta2)*(1/X.shu2-1/sh1sq)-dU(X.u);
//	Er=-E*(X.shu2-sh1sq)/ch1sq+.5*pow(X.pu,2)/(X.Delta2*ch1sq)+.5*Lzsq/(X.Delta2*ch1sq)*(1/X.shu2-1/sh1sq)+dU(X.u)/ch1sq;
	return E*(shu2-sh1sq)+Er*ch1sq-dU(u)-.5*Lzsq/(Delta*Delta)*(1/shu2-1/sh1sq);
}
double uv_orb::duturnfn(double u){//derivative of uturnfn
	double shu=sinh(u),chu=cosh(u),shu2=shu*shu,f[2];
	double Phi0=dPhiu(u,f);
	double dPhidu=Delta*(f[0]*chu*sv+f[1]*shu*cv);
	return 2*(E-Phi0)*shu*chu-(shu2+sv2)*dPhidu+Lzsq*chu/(Delta*Delta*shu2*shu);
}
double uv_orb::Rcfn(double R){//what vanishes at Rc for given Lz
	double x[2]={R,0},f[2]; dPhi(x,f);
	return Lzsq/pow(R,3)-f[0];
}
double uv_orb::Juint(double psi){
	double u=ubar-Deltau*cos(psi);
	return sqrt(MAX(0.,uturnfn(u)))*sin(psi);
}
double uv_orb::dpudE(double psi){
	double u=ubar-Deltau*cos(psi);
	return sin(psi)*pow(sinh(u),2)/sqrt(MAX(TINY,uturnfn(u)));
}
double uv_orb::dpudI3(double psi){
	double u=ubar-Deltau*cos(psi);
	return -sin(psi)/sqrt(MAX(TINY,uturnfn(u)));
}
double uv_orb::dpudLz(double psi){
	double u=ubar-Deltau*cos(psi);
	return -sin(psi)/(sqrt(MAX(TINY,uturnfn(u)))*pow(sinh(u),2));
}
double uv_orb::Phiv(double v){
	pt X(Delta,u,v);
	return Phi(X.R,X.z);
}
double uv_orb::dV(double v){//returns V(v)-V(PIH)
	return (1+shu2)*Phiv(PIH)-(shu2+pow(sin(v),2))*Phiv(v);
}
double uv_orb::vturnfn(double v){//what vanishes at apo
	double sv2=pow(sin(v),2);
	return E*sv2+Iv+dV(v)-.5*Lzsq/(Delta*Delta*sv2);
}
double uv_orb::Jvint(double psi){
	double v=PIH-(PIH-vturn)*cos(psi);
	return sqrt(MAX(0.,vturnfn(v)))*sin(psi);
}
double uv_orb::dpvdE(double psi){
	double v=PIH-(PIH-vturn)*cos(psi);
	return sin(psi)*pow(sin(v),2)/sqrt(MAX(TINY,vturnfn(v)));
}
double uv_orb::dpvdI3(double psi){
	double v=PIH-(PIH-vturn)*cos(psi);
	return sin(psi)/sqrt(MAX(TINY,vturnfn(v)));
}
double uv_orb::dpvdLz(double psi){
	double v=PIH-(PIH-vturn)*cos(psi);
	return -sin(psi)/(sqrt(MAX(TINY,vturnfn(v)))*pow(sin(v),2));
}
void uv_orb::infinity(void){
	Omegau=Omegav=Omegaphi=0; Rc1=Ju1=Jv1=1e31;
}
uv_orb::uv_orb(double Delta0,double Lz0,double Phi0,double *x,double *p){
	reset(Delta0,Lz0,Phi0,x,p);
}
void uv_orb::reset(double Delta0,double Lz0,double Phi0,double *x,double *p){
	Lz=Lz0; Lzsq=Lz*Lz; Delta=Delta0; sh1sq=pow(sinh(1),2); ch1sq=pow(cosh(1),2);
	pt X(Delta0,x,p);
	u=X.u; v=X.v; shu2=X.shu2; sv=X.sv; cv=X.cv; sv2=X.sv2;
	pv=X.pv;
	E=.5*X.p2+.5*Lzsq/pow(x[0],2)+Phi0;
	if(E>=0) infinity();//unbound particle
	else{
		sh1sq=pow(sinh(us),2); ch1sq=pow(cosh(us),2); Phiu1=Phiu(us);
		Iu=E*(X.shu2-sh1sq)-.5*pow(X.pu,2)/(X.Delta2)-.5*Lzsq/(X.Delta2)*(1/X.shu2-1/sh1sq)-dU(X.u);
		Er=-E*(X.shu2-sh1sq)/ch1sq+.5*pow(X.pu,2)/(X.Delta2*ch1sq)+.5*Lzsq/(X.Delta2*ch1sq)*(1/X.shu2-1/sh1sq)+dU(X.u)/ch1sq;
		Getu0(us);
		sh1sq=pow(sinh(umid),2); ch1sq=pow(cosh(umid),2); Phiu1=Phiu(umid);
		Iu=E*(X.shu2-sh1sq)-.5*pow(X.pu,2)/(X.Delta2)-.5*Lzsq/(X.Delta2)*(1/X.shu2-1/sh1sq)-dU(X.u);
		Er=-E*(X.shu2-sh1sq)/ch1sq+.5*pow(X.pu,2)/(X.Delta2*ch1sq)+.5*Lzsq/(X.Delta2*ch1sq)*(1/X.shu2-1/sh1sq)+dU(X.u)/ch1sq;
		Iv=.5*pow(X.pv,2)/X.Delta2-E*X.sv2+.5*Lzsq/(X.Delta2*X.sv2)-dV(X.v);
		Omegau=Omegav=Ju1=Jv1=Rc1=-1;//to signal not yet evaluated
	}
}
void uv_orb::GetTurnu(double utry){// finds apo. utry an allowed coord
	double u1,u2,fu1,fu2,fuinner,fuouter;
	double ut=uturnfn(utry),dut=duturnfn(utry),eps=.1;
	if(ut<0){
		if(dut>0){//inside peri
			uinner=utry; fuinner=ut;
			while(1){
				while(ut<0 && dut>0){
					uinner=utry; fuinner=ut; utry*=(1+eps);
					ut=uturnfn(utry); dut=duturnfn(utry);
				}
				if(ut>=0){
					u1=utry; fu1=ut; break;
				} else {
					utry/=(1+eps); eps*=.2;
					if(eps<TINY){
						//printf("too many steps from peri");
						uinner=uouter=utry; return;
					}
					utry*=(1+eps);
					ut=uturnfn(utry); dut=duturnfn(utry);
				}
			}
			while(ut>0){
				u2=utry; fu2=ut; utry*=1.1; ut=uturnfn(utry);
			}
			uouter=utry; fuouter=ut;
		}else{//outside apo
			while(1){
				while(ut<0 && dut<0){
					uouter=utry; fuouter=ut; utry/=(1+eps);
					ut=uturnfn(utry); dut=duturnfn(utry);
				}
				if(ut>=0){
					u2=utry; fu2=ut; break;
				} else {
					utry*=(1+eps); eps*=.2;
					if(eps<TINY){
						//printf("too many steps from apo");
						uinner=uouter=utry; return;
					}
					utry/=(1+eps);
					ut=uturnfn(utry); dut=duturnfn(utry);
				}
			}
			while(ut>0){
				u1=utry; fu1=ut; utry*=.9; ut=uturnfn(utry);
			}
			uinner=utry; fuinner=ut;
		}
	}else{//in allowed region
		u1=u2=utry; fu1=fu2=ut;
		while(ut>0){
			u2=utry; fu2=ut; utry*=1.1; ut=uturnfn(utry);
		}
		uouter=utry; fuouter=ut;
		ut=fu1; 
		while(ut>0){
			u1=utry; fu1=ut; utry*=.9; ut=uturnfn(utry);
		}
		uinner=utry; fuinner=ut;
	}
	uinner=Xbrent(this,&uv_orb::uturnfn,uinner,u1,fuinner,fu1,TINY,20);
	uouter=Xbrent(this,&uv_orb::uturnfn,u2,uouter,fu2,fuouter,TINY,20);
}
void uv_orb::Getu0(double utry){
	double ui,uo,dfui,dfuo,dut=duturnfn(utry);
	if(dut>0){//inside
		while(dut>0){
			ui=utry; dfui=dut; utry*=1.2; dut=duturnfn(utry);
		} uo=utry; dfuo=dut;
	}else{//outside
		while(dut<0){
			uo=utry; dfuo=dut; utry*=.8; dut=duturnfn(utry);
		} ui=utry; dfui=dut;
	}
	umid=Xbrent(this,&uv_orb::duturnfn,ui,uo,dfui,dfuo,TINY,20);
}
void uv_orb::GetRc(double R){
	double Ri,Ro,dfRi,dfRo,dfRt=Rcfn(R);
	if(dfRt>0){//inside
		while(dfRt>0){
			Ri=R; dfRi=dfRt; R*=1.2; dfRt=Rcfn(R);
		} Ro=R; dfRo=dfRt;
	}else{//outside
		while(dfRt<0){
			Ro=R; dfRo=dfRt; R*=.8; dfRt=Rcfn(R);
		} Ri=R; dfRi=dfRt;
	}
	Rc1=Xbrent(this,&uv_orb::Rcfn,Ri,Ro,dfRi,dfRo,TINY,20);
}
double uv_orb::Rc(void){
	if(Rc1==-1) GetRc(1);
	return Rc1;
}
void uv_orb::GetJu(void){//evaluates radial action
	GetTurnu(1);
	if(fabs(uouter-uinner)<TINY) Ju1=0;
	else{
		ubar=.5*(uinner+uouter);
		Deltau=.5*(uouter-uinner);
		Ju1=rtwo*Delta*Deltau*oct_int(this,&uv_orb::Juint,0,PI,4)/PI;
	}
}
double uv_orb::Ju(void){
	if(Ju1==-1) GetJu();
	return Ju1;
}
double uv_orb::GetTurnv(double vtry){// finds apo. vtry an allowed height
	double vtop=vtry,vbot=vtry;
	if(vturnfn(vtry)<0)//vtry too small; z above allowed region
		while(vturnfn(vbot)<0){
			vbot=.5*(PIH+vbot);
			if(PIH-vbot<TINY) return PIH;
		}
	else//vtry in permitted region; need to push vtop into forbidden zone
		while(vturnfn(vtop)>0){
			vtop*=.5;
			if(vtop<TINY){
				printf("vtop==TINY ");
				return 0;}
		}
	return Xbrent(this,&uv_orb::vturnfn,vbot,vtop,vturnfn(vbot),vturnfn(vtop),TINY,20);
}
void uv_orb::GetJv(void){//evaluates vert action
	vturn=GetTurnv(.5*PIH); psivmax=PIH;
	if(fabs(vturn-PIH)<TINY) Jv1=0;
	else Jv1=2*rtwo*Delta*(PIH-vturn)*oct_int(this,&uv_orb::Jvint,0,PIH,4)/PI;
}
double uv_orb::Jv(void){
	if(Jv1==-1) GetJv();
	return Jv1;
}
void uv_orb::GetFreqs(void){
	if(Ju1==-1) GetJu(); if(Jv1==-1) GetJv();
	if(vturn>PIH) vturn=PI-vturn;
	JuE=Deltau*Delta/(rtwo*PI)*oct_int(this,&uv_orb::dpudE,0,PI,4);
	JuI=Deltau*Delta/(rtwo*PI)*oct_int(this,&uv_orb::dpudI3,0,PI,4);
	JuL=Deltau*Lz/(rtwo*PI*Delta)*oct_int(this,&uv_orb::dpudLz,0,PI,4);
	/*
	 *  Update 19-02: temporary workaround to nans I was getting in Omegau-Omegav:
	 *  it happens when the apocentre is PI/2 that I get JvE,JvI,JvL=0, so Det=0.
	 *  For now I added TINY so that if vturn=PIH I can still get approximate freqs.
	 */
	JvE=rtwo*(PIH-vturn+TINY)*Delta/PI*oct_int(this,&uv_orb::dpvdE,0,PIH,4);
	JvI=rtwo*(PIH-vturn+TINY)*Delta/PI*oct_int(this,&uv_orb::dpvdI3,0,PIH,4);
	JvL=rtwo*(PIH-vturn+TINY)*Lz/(PI*Delta)*oct_int(this,&uv_orb::dpvdLz,0,PIH,4);
	double Det=JuE*JvI-JuI*JvE;
	Omegau=JvI/Det; Omegav=-JuI/Det; I3u=-JvE/Det; I3v=JuE/Det;
	Omegaphi=(JuI*JvL-JvI*JuL)/Det;
}
void uv_orb::GetThetas(double *x,double *p,double *thetau,double *thetav,double *thetaphi){
	if(Omegau==-1) GetFreqs();
	pt X(Delta,x,p);
	double psiumax=acos(MAX(-1,MIN(1,(ubar-X.u)/Deltau)));
	double JuEp=Deltau*Delta/rtwo*oct_int(this,&uv_orb::dpudE,0,psiumax,4);
	double JuIp=Deltau*Delta/rtwo*oct_int(this,&uv_orb::dpudI3,0,psiumax,4);
	double JuLp=Deltau*Lz/(rtwo*Delta)*oct_int(this,&uv_orb::dpudLz,0,psiumax,4);
	if(X.pu<0){
		JuEp=TPI*JuE-JuEp; JuIp=TPI*JuI-JuIp; JuLp=TPI*JuL-JuLp;
	}
	double psivmax = acos(MIN(1,(PIH-X.v)/(PIH-vturn)));
	double JvEp=(PIH-vturn)*Delta/rtwo*oct_int(this,&uv_orb::dpvdE,0,psivmax,4);
	double JvIp=(PIH-vturn)*Delta/rtwo*oct_int(this,&uv_orb::dpvdI3,0,psivmax,4);
	double JvLp=(PIH-vturn)*Lz/(rtwo*Delta)*oct_int(this,&uv_orb::dpvdLz,0,psivmax,4);
	if(X.pv<0){
		JvEp=TPI*JvE-JvEp; JvIp=TPI*JvI-JvIp; JvLp=TPI*JvL-JvLp;
	}
	*thetau=(JuEp+JvEp)*Omegau+(JuIp+JvIp)*I3u;
	*thetav=(JuEp+JvEp)*Omegav+(JuIp+JvIp)*I3v;
	if(*thetav<0) *thetav+=TPI; if(*thetav>TPI) *thetav-=TPI;
	*thetaphi=JuLp+JvLp;
}
#undef TINY
#undef SMALL
