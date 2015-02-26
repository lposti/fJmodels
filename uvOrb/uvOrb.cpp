/*
 * uvOrb.cpp
 *
 *  Created on: Feb 22, 2015
 *      Author: morpheus
 */

#include "uvOrb.h"
#include "Uvpt.h"
#include "Utils.h"
#include <gsl/gsl_math.h>

static const double us=58,rtwo=sqrt(2);

uvOrb::uvOrb(const double Delta0, const double Lz0, const double Phi0,
             const double *x, const double *p, Potential *pot_in){
	pot = pot_in;
	Lz=Lz0; Lzsq=Lz*Lz; Delta=Delta0; sh1sq=pow(sinh(1),2); ch1sq=pow(cosh(1),2);
	pt X(Delta0,x,p);
	u=X.u; v=X.v; shu2=X.shu2; sv=X.sv; cv=X.cv; sv2=X.sv2;
	pv=X.pv;
	E=.5*X.p2+.5*Lzsq/pow(x[0],2)+Phi0;
	if(E>=0 or isnan(E)==1) {bound=false; unbound();}//unbound particle
	else{
		bound=true;
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

void uvOrb::unbound(){
	Omegau=Omegav=Omegaphi=0; Rc1=Ju1=Jv1=1e31;
}

double uvOrb::Phiv(double v){
	pt X(Delta,u,v);
	return (*pot)(X.R,X.z);
}
double uvOrb::dV(double v){//returns V(v)-V(PIH)
	return (1+shu2)*Phiv(PIH)-(shu2+pow(sin(v),2))*Phiv(v);
}
double uvOrb::vturnfn(double v){//what vanishes at apo
	double sv2=pow(sin(v),2);
	return E*sv2+Iv+dV(v)-.5*Lzsq/(Delta*Delta*sv2);
}
double uvOrb::Jvint(double psi){
	double v=PIH-(PIH-vturn)*cos(psi);
	return sqrt(MAX(0.,vturnfn(v)))*sin(psi);
}
double uvOrb::Phiu(double u){
	pt X(Delta,u,v);
	return (*pot)(X.R,X.z);
}
double uvOrb::dU(double u){//returns U(u)-U(u0)
	return (pow(sinh(u),2)+sv2)*Phiu(u)-(sh1sq+sv2)*Phiu1;
}
double uvOrb::dPhiu(double u,double *f){//returns Phiu(u) and dPhi/dR etc
	pt X(Delta,u,v);
	double x[2]={X.R,X.z};
	f[0]=pot->dR(x[0],x[1]); f[1]=pot->dz(x[0],x[1]);
	return (*pot)(x[0],x[1]);
}
double uvOrb::uturnfn(double u){//what vanishes at peri / apo
	double shu2=pow(sinh(u),2);
	return E*(shu2-sh1sq)+Er*ch1sq-dU(u)-.5*Lzsq/(Delta*Delta)*(1/shu2-1/sh1sq);
}
double uvOrb::duturnfn(double u){//derivative of uturnfn
	double shu=sinh(u),chu=cosh(u),shu2=shu*shu,f[2];
	double Phi0=dPhiu(u,f);
	double dPhidu=Delta*(f[0]*chu*sv+f[1]*shu*cv);
	return 2*(E-Phi0)*shu*chu-(shu2+sv2)*dPhidu+Lzsq*chu/(Delta*Delta*shu2*shu);
}
double uvOrb::Juint(double psi){
	double u=ubar-Deltau*cos(psi);
	return sqrt(MAX(0.,uturnfn(u)))*sin(psi);
}

/* Wrappers for  Juint, Jvint */
double wrp_Juint(double psi, void * params){
	uvOrb * uv = (uvOrb *) params;
	return uv->Juint(psi);
}
double wrp_Jvint(double psi, void * params){
	uvOrb * uv = (uvOrb *) params;
	return uv->Jvint(psi);
}

/* Wrappers for  uturnfn, duturnfn, vturnfn */
double wrp_uturnfn(double u, void * params){
	uvOrb * uv = (uvOrb *) params;
	return uv->uturnfn(u);
}
double wrp_duturnfn(double u, void * params){
	uvOrb * uv = (uvOrb *) params;
	return uv->duturnfn(u);
}
double wrp_vturnfn(double v, void * params){
	uvOrb * uv = (uvOrb *) params;
	return uv->vturnfn(v);
}

void uvOrb::Getu0(double utry){
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

	gsl_function F;
	F.function = &wrp_duturnfn;
	F.params   = this;

	if (!(ui<uo) or isinf(duturnfn(ui))!=0 or isinf(duturnfn(uo))!=0
		or isnan(duturnfn(ui))==1 or isnan(duturnfn(uo))==1)
		    printf("Getu0: %f %f %f %f\n",ui,uo,duturnfn(ui),duturnfn(uo));
	umid = wrp_GSLroots<double>(F,ui,uo);
}

void uvOrb::GetTurnu(double utry){// finds apo. utry an allowed coord
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

	gsl_function F;
	F.function = &wrp_uturnfn;
	F.params   = this;

	if (!(uinner<u1) or !(u2<uouter)) printf("GetTurnu: %f %f %f %f\n",uinner,u1,u2,uouter);
	uinner = wrp_GSLroots<double>(F,uinner,u1);
	uouter = wrp_GSLroots<double>(F,u2,uouter);
}

double uvOrb::GetTurnv(double vtry){// finds apo. vtry an allowed height
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

	gsl_function F;
	F.function = &wrp_vturnfn;
	F.params   = this;

	if (!(vbot<vtop)) return wrp_GSLroots<double>(F,vtop,vbot);
	else return wrp_GSLroots<double>(F,vbot,vtop);
}

double uvOrb::GetJu(){//evaluates radial action
	GetTurnu(1);
	if(fabs(uouter-uinner)<TINY) Ju1=0;
	else{
		ubar=.5*(uinner+uouter);
		Deltau=.5*(uouter-uinner);

		gsl_function F;
		F.function = &wrp_Juint;
		F.params   = this;
		Ju1=rtwo*Delta*Deltau*wrp_GSLinteg<double>(F,0,PI)/PI;
	}
	return Ju1;
}

double uvOrb::Ju(){
	if(!bound) return Ju1;
	if(Ju1==-1) GetJu();
	return Ju1;
}

double uvOrb::GetJv(void){//evaluates vert action
	vturn=GetTurnv(.5*PIH); psivmax=PIH;
	if(fabs(vturn-PIH)<TINY) Jv1=0;
	else {
		gsl_function F;
		F.function = &wrp_Jvint;
		F.params   = this;
		Jv1=2*rtwo*Delta*(PIH-vturn)*wrp_GSLinteg<double>(F,0,PIH)/PI;
	}
	return Jv1;
}
double uvOrb::Jv(void){
	if(!bound) return Jv1;
	if(Jv1==-1) GetJv();
	return Jv1;
}

