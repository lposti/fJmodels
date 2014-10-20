/*
 * EtildeJ.h
 *
 *  Created on: Mar 20, 2013
 *      Author: Fermani
 */

#ifndef ETILDEJ_H_
#define ETILDEJ_H_

#include "Types.h"

class EtildeJ{
public:
	EtildeJ(double **scoeff_list,int length);
	~EtildeJ();
	int compute_surf_eigenvalues(int si,double *sc); //0 if eigenvalues are all positive, 1 otherwise
	//int closest_surf(double Jr,double Jz,double Jphi);
	//double Jr(double *coeff,double Jz,double Jphi,int *stauts);
	double Jr(int si,double Jz,double Jphi);
	vec2 Jr2(int si,double Jz,double Jphi);
	int surf_type(int i);
	double eigenvalue(int si,int i);
	double Jz(int si,double Jr,double Jphi);
	double Jphi(int si,double Jr,double Jz);
	//vec2 Jr2(double *coeff,double Jz,double Jphi,int *stauts);
	//double E(double Jr,double Jz,double Jphi,double *mJ, double *interp_coeff);
	double E(double Jr,double Jz,double Jphi,double *mJ,double *sc);
	double dEdJr(int si,double Jr,double Jz,double Jphi);
	double Omegar_Omegaz(double *gen_coeff,double Jr,double Jz,double Jphi);
	double Omegaphi_Omegaz(double *gen_coeff,double Jr,double Jz,double Jphi);
	double Omegar_Omegaphi(double *gen_coeff,double Jr,double Jz,double Jphi);
	double Omegaz_Omegaphi(double *gen_coeff,double Jr,double Jz,double Jphi);
	double Omegaz_Omegar(double *gen_coeff,double Jr,double Jz,double Jphi);
	double Omegaphi_Omegar(double *gen_coeff,double Jr,double Jz,double Jphi);
	void testVSIso();
private:
	int lib_size_;
	double **surf_library_;
	int *surf_type_; //0 for ellipsoid, 1 for hyperboloid
	double **surf_eigenval_;
};


#endif /* ETILDEJ_H_ */
