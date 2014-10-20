/*
 * numericPot.h
 *
 *  Created on: 13/mag/2014
 *      Author: morpheus
 */

#ifndef NUMERICPOT_H_
#define NUMERICPOT_H_

#include "potential_obj.h"

class numericPot: public potential_obj {
public:
	numericPot(int npar,double *coeff,char *fname);
	virtual ~numericPot();

	void readPotential();
	void getCurrentPot();
	double phi(double R,double z);
	double dphidR(double R,double z);
	double dphidz(double R,double z);
	double PhiR(double R);
	double dPhiR(double R);
	double dPhi(double* x2,double *f2);

private:
	int nr_, npoly_, ngauss_;
	double *ar_, **rhl_, **phil_, **Pr_, **Pr2_;

	int potleg_();
	double evpot_(double R2,double z);
	double evfor_(double R2,double z);
	void intpo2_(double r,double *phip);
	void intfor_(double r,double *dphip);
	double dPhir_(double,double);
	double dPhitheta_(double,double);
	void legend_(double *,double,int);
	double dlegend_(double,int);
	template<typename T> void evenlegend_(T *pol,T c);
};

/*
 * template for evenlegend
 */
template<typename T> void numericPot::evenlegend_(T *pol,T c){
// evaluates even legendre polys up to l=2*(npoly-1) at c -------------
	T c2=c*c;
	pol[0]=1; if(npoly_<2) return;
	pol[1]=1.5*c2-.5;
	for(int np=2; np<npoly_; np++){
		int l=2*(np-1), l2=2*l;
		pol[np]=-pol[np-2]*l*(l-1)/(T)((l2+1)*(l2-1))+
			  pol[np-1]*(c2-(l2*l+l2-1)/(T)((l2-1)*(l2+3)));
		pol[np]*=(l2+1)*(l2+3)/(T)((l+1)*(l+2));
	}
}

#endif /* NUMERICPOT_H_ */
