/*
 * EtildeJ.cpp
 *
 *  Created on: Mar 20, 2013
 *      Author: Fermani
 */
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "EtildeJ.h"
#include "Iso_check.h"
#include "Types.h"
#include "Structure.h"
#include "tools.h"

using namespace std;

using std::cout;
using std::cerr;
using std::string;


//double coeff[2]={1.0e12, 1.0};
//Iso_check *tiso;


EtildeJ::EtildeJ(double **scoeff_list,int length){
	surf_library_=scoeff_list;
	lib_size_=length;
	surf_type_=new int[length];
	surf_eigenval_=new double* [length];
	for(int j=0;j<length;j++){
		surf_eigenval_[j]=new double[3];
		surf_type_[j]=compute_surf_eigenvalues(j,scoeff_list[j]);
	}
}

EtildeJ::~EtildeJ(){}

int EtildeJ::compute_surf_eigenvalues(int si,double *sc){
	int i;
	double **M=new double* [3];
	for(i=0;i<3;i++) M[i]=new double[3];

	M[0][0]=sc[4]; M[0][1]=0.5*sc[7]; M[0][2]=0.5*sc[9];
	M[1][0]=0.5*sc[7]; M[1][1]=sc[5]; M[1][2]=0.5*sc[8];
	M[2][0]=0.5*sc[9]; M[2][1]=0.5*sc[8]; M[2][2]=sc[6];

	double eig1,eig2,eig3;
	eigenvalues_33(&eig1,&eig2,&eig3,M);

	surf_eigenval_[si][0]=eig1; surf_eigenval_[si][1]=eig2; surf_eigenval_[si][2]=eig3;
	for(i=0;i<3;i++){ delete[] M[i];} delete[] M;
	if(eig1>=0 && eig2>=0 && eig3>=0) return 0;
	else return 1;
}

int EtildeJ::surf_type(int i){ return surf_type_[i];}

double EtildeJ::eigenvalue(int si,int i){ return surf_eigenval_[si][i];}

double EtildeJ::E(double Jr,double Jz,double Jphi,double *mJ,double *sc){
	double approx_value, wl, wu;
	Jphi=fabs(Jphi);
	int i,l,u;
	double Jrl,Jru;


	// store Jr values in an array
	double **Jr_array;
	Jr_array= new double* [lib_size_];
	for(i=0;i<lib_size_;i++) Jr_array[i]=new double[2];
	for(i=0;i<lib_size_;i++){
		Jr_array[i][0]=surf_library_[i][0];
		Jr_array[i][1]=this->Jr(i,Jz,Jphi);
	}
	// Sort Jr_array: aimed order is from J_max down to J_min
	int ite=0,iteint=0,iteint_temp;
	double t_maxJr;
	int *Jrorder,eff_length;
	Jrorder=new int[lib_size_];

	double **arrayres,**arraytemp;
	arrayres= new double* [lib_size_]; arraytemp= new double* [lib_size_];
	for(i=0;i<lib_size_;i++){ arrayres[i]=new double[2]; arraytemp[i]=new double[2];}

	ite=0;
	Jrorder[ite]=0; t_maxJr=Jr_array[0][1];
	for(i=1;i<lib_size_;i++) if(Jr_array[i][1]<2000.0 && Jr_array[i][1]>t_maxJr){ Jrorder[ite]=i; t_maxJr=Jr_array[i][1];}
	t_maxJr=Jr_array[Jrorder[ite]][1];
	ite++;

	iteint=0;
	for(i=0;i<lib_size_;i++) if(Jr_array[i][1]>=0 && Jr_array[i][1]<t_maxJr && Jr_array[i][1]>0){
		arraytemp[iteint][0]=i;
		arraytemp[iteint][1]=Jr_array[i][1];
		iteint++;
	}

	while(iteint>1){
		//cout << ite << " " << iteint << "\n";
		Jrorder[ite]=arraytemp[0][0];
		t_maxJr=arraytemp[0][1];
		for(i=1;i<iteint;i++) if(arraytemp[i][1]>t_maxJr){ Jrorder[ite]=arraytemp[i][0]; t_maxJr=arraytemp[i][1];}
		t_maxJr=Jr_array[Jrorder[ite]][1];
		ite++;

		iteint_temp=iteint;
		iteint=0;
		for(i=0;i<iteint_temp;i++) if(arraytemp[i][1]<t_maxJr){
			arraytemp[iteint][0]=arraytemp[i][0];
			arraytemp[iteint][1]=arraytemp[i][1];
			iteint++;
		}
	}
	Jrorder[ite]=arraytemp[0][0];
	ite++;
	eff_length=ite;

    // find closest value to Jr in the Jr_array
	double jrdiscp,jrdiscp_pp;
	int ind=0;
	if(Jr>Jr_array[Jrorder[0]][1] || Jr<Jr_array[Jrorder[eff_length-1]][1]){ //If Jr outside array range, can't interpolate
		//cout << "I cannot perform interpolation: Jr outside table. Take table edge values instead. Jr=" << Jr << ", Jr_edges=" << Jr_array[Jrorder[0]][1] << " " << Jr_array[Jrorder[eff_length-1]][1] << ".\n";
		if(Jr>Jr_array[Jrorder[0]][1]){
			l=Jrorder[0]; u=Jrorder[0]; wl=1.0; wu=0.0;
		}
		else{ l=Jrorder[eff_length-1]; u=Jrorder[eff_length-1]; wl=0.0; wu=1.0;}
	}
	else{
		jrdiscp=fabs(Jr_array[Jrorder[ind]][1]-Jr);
		jrdiscp_pp=fabs(Jr_array[Jrorder[ind+1]][1]-Jr);
		while(jrdiscp_pp<jrdiscp && ind<eff_length-2){
			ind++;
			jrdiscp=fabs(Jr_array[Jrorder[ind]][1]-Jr);
			jrdiscp_pp=fabs(Jr_array[Jrorder[ind+1]][1]-Jr);
		}
		// determine the enclosing array values
		if(Jr_array[Jrorder[ind]][1]-Jr>0){ Jru=Jr_array[Jrorder[ind]][1]; Jrl=Jr_array[Jrorder[ind+1]][1]; u=Jrorder[ind]; l=Jrorder[ind+1];}
		else{ Jru=Jr_array[Jrorder[ind-1]][1]; Jrl=Jr_array[Jrorder[ind]][1]; u=Jrorder[ind-1]; l=Jrorder[ind];}
		wl=(Jr-Jrl)/(Jru-Jrl); wu=(Jru-Jr)/(Jru-Jrl);
		if(wl>1){// || wu >1){
			cerr << "Shout for help! [weights bigger than 1]\n"; cout << Jr << " " << Jrl << " " << Jru << " " << wl << " " << wu << "\n";
			exit(1);
		}
	}
	approx_value=wl*Jr_array[l][0]+wu*Jr_array[u][0];
	//tiso=new Iso_check(2,&coeff[0]);
	//cout << "[E check] "

	double S_a,S_b,S_c,q,S_b24ac, mjt1,mjt2,mjl,mju;

	sc[0]=approx_value;
	for(i=1;i<10;i++) sc[i]=wl*surf_library_[l][i]+wu*surf_library_[u][i];

	S_a=sc[4]+sc[5]+sc[6]+sc[7]+sc[8]+sc[9];
	S_b=sc[1]+sc[2]+sc[3];
	S_c=-sc[0];
	if(pow(S_b,2.0)-4.0*S_a*S_c>=0){
		q=-0.5*(S_b+sign(S_b)*sqrt(pow(S_b,2.0)-4.0*S_a*S_c));
		mjt1=q/S_a; mjt2=S_c/q;
		if(mjt1<mjt2 && mjt1>=0) *mJ=mjt1;
		else *mJ=mjt2;
	}
	else{
		*mJ=-S_b/(2.0*S_a);
	}

	for(i=0;i!=lib_size_;i++){
		delete[] Jr_array[i]; delete[] arrayres[i]; delete[] arraytemp[i];
	}
	delete[] Jr_array; delete[] arrayres; delete[] arraytemp; delete[] Jrorder;

	return approx_value;
}


double EtildeJ::Jr(int si,double Jz,double Jphi){
	double A,B,C;
	double q,J01,J02,linear_sol;

	A=surf_library_[si][4]; B=surf_library_[si][1]+surf_library_[si][7]*Jz+surf_library_[si][9]*Jphi;
	C=surf_library_[si][2]*Jz+surf_library_[si][3]*Jphi+surf_library_[si][5]*pow(Jz,2.0)+surf_library_[si][6]*pow(Jphi,2.0)+surf_library_[si][8]*Jz*Jphi-surf_library_[si][0];

	linear_sol=(surf_library_[si][0]-surf_library_[si][2]*Jz-surf_library_[si][3]*Jphi)/surf_library_[si][1];


	if(pow(B,2.0)-4.0*A*C>=0){
		q=-0.5*(B+sign(B)*sqrt(pow(B,2.0)-4.0*A*C));
		J01=q/A; J02=C/q;

		if(fabs(this->eigenvalue(si,1))<fabs(this->eigenvalue(si,2))) return min(J01,J02);
		else return max(J01,J02);

	}
	else{
		return -B/(2.0*A);
	}

}

vec2 EtildeJ::Jr2(int si,double Jz,double Jphi){
	double A,B,C;
	double q,J01,J02,linear_sol;
	vec2 res;

	A=surf_library_[si][4]; B=surf_library_[si][1]+surf_library_[si][7]*Jz+surf_library_[si][9]*Jphi;
	C=surf_library_[si][2]*Jz+surf_library_[si][3]*Jphi+surf_library_[si][5]*pow(Jz,2.0)+surf_library_[si][6]*pow(Jphi,2.0)+surf_library_[si][8]*Jz*Jphi-surf_library_[si][0];

	linear_sol=(surf_library_[si][0]-surf_library_[si][2]*Jz-surf_library_[si][3]*Jphi)/surf_library_[si][1];


	if(pow(B,2.0)-4.0*A*C>=0){
		q=-0.5*(B+sign(B)*sqrt(pow(B,2.0)-4.0*A*C));
		J01=q/A; J02=C/q;
		if(fabs(this->eigenvalue(si,1))<fabs(this->eigenvalue(si,2))){ res[0]=min(J01,J02); res[1]=max(J01,J02);}
		else{ res[0]=max(J01,J02); res[1]=min(J01,J02);}

	}
	else{
		res[0]=-B/(2.0*A);
		res[1]=-B/(2.0*A);
	}
	return res;
}

double EtildeJ::Jz(int si,double Jr,double Jphi){
	double A,B,C,Jz=9999.9;
	double q,J01,J02;

	A=surf_library_[si][5]; B=surf_library_[si][2]+surf_library_[si][7]*Jr+surf_library_[si][8]*Jphi;
	C=surf_library_[si][1]*Jr+surf_library_[si][3]*Jphi+surf_library_[si][4]*pow(Jr,2.0)+surf_library_[si][6]*pow(Jphi,2.0)+surf_library_[si][9]*Jr*Jphi-surf_library_[si][0];

	if(pow(B,2.0)-4.0*A*C>=0){
		q=-0.5*(B+sign(B)*sqrt(pow(B,2.0)-4.0*A*C));
		J01=q/A; J02=C/q;
		if(J01<J02) Jz=J01;
		else Jz=J02;
	}
	return Jz;
}

double EtildeJ::Jphi(int si,double Jr,double Jz){
	double A,B,C,Jphi=9999.9;
	double q,J01,J02;

	A=surf_library_[si][6]; B=surf_library_[si][3]+surf_library_[si][8]*Jz+surf_library_[si][9]*Jr;
	C=surf_library_[si][1]*Jr+surf_library_[si][2]*Jz+surf_library_[si][4]*pow(Jr,2.0)+surf_library_[si][5]*pow(Jz,2.0)+surf_library_[si][7]*Jr*Jz-surf_library_[si][0];

	if(pow(B,2.0)-4.0*A*C>=0){
		q=-0.5*(B+sign(B)*sqrt(pow(B,2.0)-4.0*A*C));
		J01=q/A; J02=C/q;
		if(J01<J02) Jphi=J01;
		else Jphi=J02;
	}
	return Jphi;
}

double EtildeJ::dEdJr(int si,double Jr,double Jz,double Jphi){
	return surf_library_[si][1]+2.0*surf_library_[si][4]*Jr+surf_library_[si][7]*Jz+surf_library_[si][9]*Jphi;
}

double EtildeJ::Omegar_Omegaz(double *gen_coeff,double Jr,double Jz,double Jphi){
	return (gen_coeff[1]+2.0*gen_coeff[4]*Jr+gen_coeff[7]*Jz+gen_coeff[9]*Jphi)/(gen_coeff[2]+2.0*gen_coeff[5]*Jz+gen_coeff[7]*Jr+gen_coeff[8]*Jphi);
}

double EtildeJ::Omegaphi_Omegaz(double *gen_coeff,double Jr,double Jz,double Jphi){
	return (gen_coeff[3]+2.0*gen_coeff[6]*Jphi+gen_coeff[8]*Jz+gen_coeff[9]*Jr)/(gen_coeff[2]+2.0*gen_coeff[5]*Jz+gen_coeff[7]*Jr+gen_coeff[8]*Jphi);
}

double EtildeJ::Omegar_Omegaphi(double *gen_coeff,double Jr,double Jz,double Jphi){
	return (gen_coeff[1]+2.0*gen_coeff[4]*Jr+gen_coeff[7]*Jz+gen_coeff[9]*Jphi)/(gen_coeff[3]+2.0*gen_coeff[6]*Jphi+gen_coeff[8]*Jz+gen_coeff[9]*Jr);
}

double EtildeJ::Omegaz_Omegaphi(double *gen_coeff,double Jr,double Jz,double Jphi){
	return 1.0/this->Omegaphi_Omegaz(gen_coeff,Jr,Jz,Jphi);
}

double EtildeJ::Omegaz_Omegar(double *gen_coeff,double Jr,double Jz,double Jphi){
	return 1.0/this->Omegar_Omegaz(gen_coeff,Jr,Jz,Jphi);
}

double EtildeJ::Omegaphi_Omegar(double *gen_coeff,double Jr,double Jz,double Jphi){
	return 1.0/this->Omegar_Omegaphi(gen_coeff,Jr,Jz,Jphi);
}

void EtildeJ::testVSIso(){
	string input("test_ham_approx.dat");
	const char *in_p=(const char*)input.c_str();
	ofstream toT;

	double Iso_par[2]={1.0e12,1.0};
	Iso_check *tiso;
	tiso = new Iso_check(2,Iso_par);

	int grid=10,pt_exc=0;
	double sc[10];
	double J0,J1,J2;
	double discrepancy=0,sigma_d=0;
	double e_fit,mJ,maxJ;
	maxJ=tiso->L(-4.8,0.); // in Pauls' unit this corresponds to E=-0.05

	//toT.open(in_p);
	for(int i=0;i!=int(grid);i++){
		for(int j=0;j!=int(grid);j++){
			for(int k=0;k!=int(grid);k++){
				J0=double(i)*(maxJ/double(grid))+1.0;
				J1=double(j)*(maxJ/double(grid))+1.0;
				J2=double(k)*(maxJ/double(grid))+1.0;
				//cout << "[Jr,Jz,Jphi] " << J0 << " " << J1 << " " << J2 << "\n";
				e_fit=this->E(J0,J1,J2,&mJ,&sc[0]);
				if(e_fit<0 && tiso->E(J0,J1,fabs(J2))<=-0.05*unit2){
					if(J2!=0 && (e_fit+2.0)/2.0==0.5*e_fit+1.0){
						//toT << J0 << " " << J1 << " " << J2 << " " << tiso->E(J0,J1,fabs(J2)) << " " << e_fit << " " << (tiso->E(J0,J1,fabs(J2))-e_fit)/tiso->E(0.,0.,0.) << "\n";
						discrepancy+=(tiso->E(J0,J1,fabs(J2))-e_fit)/tiso->E(0.,0.,0.);
						sigma_d+=pow((tiso->E(J0,J1,fabs(J2))-e_fit)/tiso->E(0.,0.,0.),2.0);
					}
					else pt_exc++;
				}
				else pt_exc++;
			}
		}
	}
	//cout << "[grid,pt_exc,maxJ] " << grid << " " << pt_exc << " " << maxJ << "\n";
	//cout << "[disc, sigma_d] " << discrepancy << " " << sigma_d << "\n";

	grid=grid*grid*grid-pt_exc;
	discrepancy=discrepancy/double(grid);
	sigma_d=sqrt(fabs(sigma_d/double(grid)-pow(discrepancy,2.0)))/sqrt(double(grid));

	//cout << "[grid,pt_exc,maxJ] " << grid << " " << pt_exc << " " << maxJ << "\n";
	//cout << "[disc, sigma_d] " << discrepancy << " " << sigma_d << "\n";

	cout << "Discrepancy in E (linear interp. " << lib_size_ << " surfaces): " << discrepancy << " +- " << sigma_d << "\n";
	//toT.close();
}

