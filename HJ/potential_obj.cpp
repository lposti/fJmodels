/*
 * potential.cpp
 *
 *  Created on: Mar 15, 2013
 *      Author: Fermani
 */
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "Structure.h"
#include "potential_obj.h"
#include "falPot.h"
#include "Units.h"
#include "isopotF.h"

using std::cout;
using std::cin;

double Rzcutoff=rsup;

potential_obj::potential_obj(int pot_index,int npar,double *coeff,char *fname){
	ifstream file;

	whichpotential_=pot_index;
	npar_=npar;
	fname_=fname;
	if(whichpotential_==1){
		file.open(fname_);
		PhiWD_=new GalaxyPotential(file);
		file.close();
	}
	for(int i=0;i<npar;i++) a_[i]=coeff[i];
	if(npar>2 && whichpotential_==4) isopot_init(coeff[1],coeff[2]); // should be b and q
}

potential_obj::~potential_obj(){
	delete &whichpotential_;
	delete &npar_;
	delete &a_[0];
}

double potential_obj::phi(double R,double z){
    double Hern_pot(double R,double z,double *coeff);
    double Iso(double R,double z,double *coeff);
    double MN(double R,double z,double *coeff);
    double logPot(double R,double z,double *coeff);
    double myIso(double R,double z,double *coeff);
    double res;

    double coeff[npar_];
    for(int i=0;i<npar_;i++) coeff[i]=a_[i];

    switch(whichpotential_){
    	case 0:
    		res=Hern_pot(R,z,&coeff[0]);
    		break;
    	case 1:
    		res=unit2*(*PhiWD_)(R,z);
    		break;
    	case 2:
    		res=MN(R,z,&coeff[0]);
    		break;
    	case 3:
    		res=Iso(R,z,&coeff[0]);
    		break;
    	case 4:
    		res=unit2*U.G*coeff[0]*isopot(R,z);
    		break;
    	/* 12/05/14
    	 * LP edit: added isothermal sphere
    	 */
    	case 5:
    		res=myIso(R,z,&coeff[0]);//logPot(R,z,&coeff[0]);
    		break;
    }
    return res;
}

double potential_obj::dphidR(double R,double z){
	double res, a1, b, M, R2azb2, temp, f2[2], rho0, sig, rmax;
	switch(whichpotential_){
		case 0:
			M=a_[0]; a1=a_[1];
			rho0=M*pow(a1+Rzcutoff,2.0)/(2.0*a1*a1*a1*Pi*Rzcutoff*Rzcutoff);
			if(z==0 && R==0) res= 2.0*a1*U.G*Pi*rho0;
			else{
				if(R<SMALL) res=0.0;
				else res= 2.0*a1*U.G*Pi*R*rho0/(sqrt(R*R+z*z)*pow(1.0+sqrt(R*R+z*z)/a1,2.0));
			}

			if(R>Rzcutoff || z>Rzcutoff){
				res=0.0;
			}
			break;
		case 1:
			temp=unit2*(*PhiWD_)(R,z,f2[0],f2[1]);
			f2[0]*=unit2; f2[1]*=unit2;
			res=f2[0];
			break;
		case 2:
			M=a_[0], a1=a_[1], b=a_[2];
			R2azb2= pow(R,2.0)+pow(a1+sqrt(z*z+b*b),2.0);
			if(R==0) res=0.0;
			else res=U.G*M*R/pow(R2azb2,1.5);
			break;
		case 3:
			M=a_[0]; b=a_[1];
			res = U.G*M*R/(sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0))*pow(b+sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0)),2.0));
			break;
		case 4:
			M=a_[0];
			res=-isoforceFR(R,z);
			res*=U.G*M;
			break;
		/* 12/05/14
		 * LP edit: added isothermal sphere
		 */
	   	case 5:
	   		M=a_[0]; b=a_[1];
	   		//rmax=75.; sig = sqrt(6.67e-8*M/(2.*rmax*3.08e21))/1e5;
	   		//res = 2*R/(pow(b,2.0)+pow(R,2.0)+pow(z,2.0));//(sqrt(pow(R,2.0)+pow(z,2.0))*(b+sqrt(pow(R,2.0)+pow(z,2.0))));

	   		res = 4.5*M*R/(sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0))*pow(b+sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0)),2.0));
	   		break;
	}
	return unit2*res;
}
double potential_obj::dphidz(double R,double z){
	double res, a1, b, M, R2azb2, temp, f2[2], rho0, sig, rmax;
	switch(whichpotential_){
		case 0:
			M=a_[0]; a1=a_[1];
			rho0=M*pow(a1+Rzcutoff,2.0)/(2.0*a1*a1*a1*Pi*Rzcutoff*Rzcutoff);
			if(z==0 && R==0) res= 2.0*a1*U.G*Pi*rho0;
			else{
				if(z<SMALL) res=0.0;
				else res= 2.0*a1*U.G*Pi*z*rho0/(sqrt(R*R+z*z)*pow(1.0+sqrt(R*R+z*z)/a1,2.0));
			}
			if(R>Rzcutoff || z>Rzcutoff){
				res=0.0;
			}
			break;
		case 1:
			temp=unit2*(*PhiWD_)(R,z,f2[0],f2[1]);
			f2[0]*=unit2; f2[1]*=unit2;
			res=f2[1];
			break;
		case 2:
			M=a_[0], a1=a_[1], b=a_[2];
			R2azb2= pow(R,2.0)+pow(a1+sqrt(z*z+b*b),2.0);
			if(z==0) res=0.0;
			else res=U.G*M*z*(a1+sqrt(b*b+z*z))/(sqrt(b*b+z*z)*pow(R2azb2,1.5));
			break;
		case 3:
			M=a_[0]; b=a_[1];
			res = U.G*M*z/(sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0))*pow(b+sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0)),2.0));
			break;
		case 4:
			M=a_[0];
			res=-isoforceFz(R,z);
			res*=U.G*M;
			break;
		/* 12/05/14
		 * LP edit: added isothermal sphere
		 */
	 	case 5:
	 		M=a_[0]; b=a_[1];
	 		//rmax=75.; sig = sqrt(6.67e-8*M/(2.*rmax*3.08e21))/1e5;
	 		//res = 2*z/(pow(b,2.0)+pow(R,2.0)+pow(z,2.0));//(sqrt(pow(R,2.0)+pow(z,2.0))*(b+sqrt(pow(R,2.0)+pow(z,2.0))));

	 		res = 4.5*M*z/(sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0))*pow(b+sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0)),2.0));
	   		break;
	}

	return unit2*res;
}

double potential_obj::PhiR(double R){ return this->phi(R,0.);}

double potential_obj::dPhiR(double R){return this->dphidR(R,0.);}

double potential_obj::dPhi(double *x2,double *f2){
	f2[0]=this->dphidR(x2[0],x2[1]);
	f2[1]=this->dphidz(x2[0],x2[1]);
	return this->phi(x2[0],x2[1]);
}

double potential_obj::M_coeff(){return a_[0];}

double potential_obj::b_coeff(){return a_[1];}

// Hernquist Potential:

double Hern_pot(double R,double z,double *coeff){
	double M=coeff[0], a=coeff[1], res, max_r=Rzcutoff*sqrt(2.);
	double rho0=M*pow(a+max_r,2.0)/(2.0*a*a*a*Pi*max_r*max_r);

    res= -2.0*Pi*U.G*rho0*pow(a,2.0)/(1.0+sqrt(R*R+z*z)/a);

	if(R>Rzcutoff || z>Rzcutoff){
		res=0.0;
	}
	return unit2*res;
}

// Miyamoto-Nagai Potential:

double MN(double R,double z,double *coeff){
	double M=coeff[0], a=coeff[1], b=coeff[2];
	double R2azb2= pow(R,2.0)+pow(a+sqrt(z*z+b*b),2.0),res;
	res= -U.G*M/sqrt(R2azb2);

	return unit2*res;
}

// Isochrone Potential:

double Iso(double R,double z,double *coeff){
	double M=coeff[0], b=coeff[1];
	double res = -U.G*M/(b+sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0)));
	return unit2*res;
}

/*
 * 12/05/14
 * LP edit: logarithmic potential
 */
double myIso(double R,double z,double *coeff){
	double M=coeff[0], b=coeff[1];
	double res = -M/(b+sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0)));
	return 430.*res;
}
double logPot(double R,double z,double *coeff){
	double M=coeff[0], b=coeff[1];
	double rmax=75.; //sig = sqrt(6.67e-8*M*1.989e33/(2.*rmax*3.08e21))/1e5;
	double res = 2.*(log(sqrt(pow(b,2.0)+pow(R,2.0)+pow(z,2.0)))-log(rmax));
	return res;
}
