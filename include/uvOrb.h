/*
 * uvOrb.h
 *
 *  Created on: Feb 22, 2015
 *      Author: L. Posti
 */

#ifndef INCLUDE_UVORB_H_
#define INCLUDE_UVORB_H_

#include "Potential.h"

class uvOrb{
public:
	uvOrb(const double,const double,const double,const double*,const double*,Potential *);

	Potential *pot;
	double Delta,u,v,pv,Iu,Iv;
	double E,Er,Lz,Omegau,Omegav,Omegaphi;
	double uinner,uouter,umid,vturn;

	double Rc(void);
	double Ju(void);
	double Jv(void);
	void GetThetas(double*,double*,double*,double*,double*);
	void GetFreqs(void);
	double Phiu(double);
	double dU(double);
	double dPhiu(double,double*);
	double uturnfn(double);
	double duturnfn(double);
	double Rcfn(double);
	double Juint(double);
	double dpudE(double);
	double dpudI3(double);
	double dpudLz(double);
	double Phiv(double);
	double dV(double);
	double vturnfn(double);
	double Jvint(double);
	double dpvdE(double);
	double dpvdI3(double);
	double dpvdLz(double);

private:
	bool bound;
	double Rc1,Ju1,Jv1,I3u,I3v;
	double JuE,JuI,JuL,JvE,JvI,JvL;
	double Lzsq,psivmax,Deltau;
	double shu2,sh1sq,ch1sq,Phiu1,ubar;
	double sv,cv,sv2;
	void GetTurnu(double);
	double GetJu(void);
	double GetTurnv(double);
	double GetJv(void);
	void Getu0(double);
	void GetRc(double);
	void unbound();
};



#endif /* INCLUDE_UVORB_H_ */
