/*
 * Iso_check.cpp
 *
 *  Created on: Mar 16, 2013
 *      Author: Fermani
 */

#include <stdlib.h>
#include <cmath>
#include "Iso_check.h"
#include "Structure.h"
#include "Units.h"

using namespace std;

using std::cout;

Iso_check::Iso_check(int npar,double *coeff){
	npar_=npar;
	for(int i=0;i<npar;i++) a_[i]=coeff[i];
}

Iso_check::~Iso_check(){}

double Iso_check::Jr(double E,double Jz,double Jphi){
	E/=unit2; Jz/=unit; Jphi/=unit;
	double M=a_[0], b=a_[1];
	double L=Jz+fabs(Jphi);
	cout << "   E: " << E << "   L: " << L << "  JR= " << U.G*M/sqrt(-2.0*E) << " " << -0.5*(L+sqrt(pow(L,2.0)+4.0*U.G*M*b)) << endl;
	return unit*(U.G*M/sqrt(-2.0*E)-0.5*(L+sqrt(pow(L,2.0)+4.0*U.G*M*b)));
}

double Iso_check::Jphi(double E,double Jr,double Jz){
	E/=unit2; Jr/=unit; Jz/=unit;
	double M=a_[0], b=a_[1];
	return unit*(-Jz+1.0/(E*(E*Jr*Jr+0.5*U.G*U.G*M*M))*(-0.353553*pow(U.G,4.0)*sqrt(-E/pow(U.G*M,2.0))*pow(M,4.0)+E*E*(-pow(Jr,3.0)+b*U.G*Jr*M)+pow(U.G*M,2.0)*E*(-0.5*Jr-0.707107*Jr*Jr*sqrt(-E/pow(U.G*M,2.0))-0.707107*b*U.G*sqrt(-E/pow(U.G*M,2.0))*M)));
}

double Iso_check::E(double Jr,double Jz,double Jphi){
	Jr/=unit; Jz/=unit; Jphi/=unit;
	double M=a_[0], b=a_[1];
	double L=Jz+fabs(Jphi);
	return unit2*(-0.5*pow(U.G*M,2.0)/pow(Jr+0.5*(L+sqrt(L*L+4.0*U.G*M*b)),2.0));
}

double Iso_check::L(double E,double Jr){
	E/=unit2; Jr/=unit;
	double M=a_[0], b=a_[1], gmb=U.G*M*b, jr2E=Jr*sqrt(-2.0*E);
	return unit*(gmb*sqrt(-2.0*E)/(jr2E-U.G*M)-(jr2E-U.G*M)/sqrt(-2.0*E));
}

double Iso_check::DeltaJr(double **scoeff,int sindex,double Jz,double Jphi){
	double Jr_fit(double **scoeff,int sindex,double Jz,double Jphi,int *status);

	int status;
	double Jr_true=this->Jr(scoeff[sindex][0],Jz,Jphi);
	double Jr_quad=Jr_fit(scoeff,sindex,Jz,Jphi,&status);

	//cout << "[Jz,Jphi] " << Jz << " " << Jphi << "\n";
	//cout << "[Jr] " << Jr_true << " " << Jr_quad << " " << status << "\n";
	if(status==0) return (Jr_true-Jr_quad)/(Jr_true+Jz+fabs(Jphi));
	else return 0.;
}

double Iso_check::DeltaE(double **scoeff,int sindex,double Jz,double Jphi){
	double E_quad(double **scoeff,int sindex,double Jr,double Jz,double Jphi);
	double Jr=this->Jr(scoeff[sindex][0],Jz,Jphi);
	//cout << "[E] " << scoeff[sindex][0] << " " << E_quad(scoeff,sindex,Jr,Jz,Jphi) << "\n";
	return (scoeff[sindex][0]-E_quad(scoeff,sindex,Jr,Jz,Jphi))/scoeff[sindex][0];
}

double Iso_check::rho(double r){
	double M=a_[0],b=a_[1],a=sqrt(b*b+r*r);
	return M*(3.*(b+a)*a*a-r*r*(b+3.*a))/(4.0*Pi*pow(b+a,3.0)*pow(a,3.0));
}

double Jr_fit(double **scoeff,int sindex,double Jz,double Jphi,int *status){
	double A,B,C,Jr=9999.9;
	double q,J01,J02;

	A=scoeff[sindex][4]; B=scoeff[sindex][1]+scoeff[sindex][7]*Jz+scoeff[sindex][9]*Jphi;
	C=scoeff[sindex][2]*Jz+scoeff[sindex][3]*Jphi+scoeff[sindex][5]*pow(Jz,2.0)+scoeff[sindex][6]*pow(Jphi,2.0)+scoeff[sindex][8]*Jz*Jphi-scoeff[sindex][0];

	if(pow(B,2.0)-4.0*A*C>=0){
		q=-0.5*(B+sign(B)*sqrt(pow(B,2.0)-4.0*A*C));
		J01=q/A; J02=C/q;
		if(J01<J02) Jr=J01;
		else Jr=J02;
		*status=0;
	}
	else *status=1;
	return Jr;
}

double E_quad(double **scoeff,int sindex,double Jr,double Jz,double Jphi){
	return scoeff[sindex][1]*Jr+scoeff[sindex][2]*Jz+scoeff[sindex][3]*Jphi+scoeff[sindex][4]*pow(Jr,2.0)+scoeff[sindex][5]*pow(Jz,2.0)+scoeff[sindex][6]*pow(Jphi,2.0)+scoeff[sindex][7]*Jr*Jz+scoeff[sindex][8]*Jz*Jphi+scoeff[sindex][9]*Jr*Jphi;
}
