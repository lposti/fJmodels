/*
 * ini_potLEG.h
 *
 *  Created on: 05/mar/2014
 *      Author: morpheus
 */

#ifndef INI_POTLEG_H_
#define INI_POTLEG_H_

void intpo2_ini(double r,double *phip,int npoly);
void intpo2_ini(double r,double *phip,double *dphip,double *d2phip,int npoly);
double dPhir_ini(double r,double theta,double *d2Phidr2);
double dPhitheta_ini(double r, double c);
double d2Phitheta_ini(double r, double c);
double d2PhiR_ini(double R);
double d2Phiz_ini(double R);
double dPhiR_ini(double R, double z);
double dPhiz_ini(double R, double z);
double dPhi_ini(double *x2,double *f2);
void getfreqs_ini(double R,double *kappa,double *nu,double *Omega);


#endif /* INI_POTLEG_H_ */
