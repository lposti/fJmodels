/*
 * isodf4.h
 *
 *  Created on: 05/feb/2014
 *      Author: morpheus
 */

#ifndef ISODF4_H_
#define ISODF4_H_

#include "HJ/surf_lib.h"
#include "HJ/EtildeJ.h"

//#define ISOCHRONE
//#define NFW
#define HERNQUIST
//#define ISOTHERMAL

extern double *EofJ,*LzofJ;
extern int rep;

extern surf_lib *sl;
extern double **surf_LIB;
extern EtildeJ *EtJ;

void setdf(double aphi_in,double az_in);
double df_iso(double Et);
void build_df_grid(void);
void build_EofJ(void);
void delete_df_grid(void);
double df_I(double Jr,double Lz,double Jz,double R);
double df_H(double Jr,double Lz,double Jz,double *x, double *v);
double expdf_I(double j,double Jr,double Lz,double Jz,double R,
             double deltar,double deltaphi,double deltaz);
double expdf_H(double j,double Jr,double Lz,double Jz,
			 double *x, double *v,
             double deltar,double deltaphi,double deltaz);
double df(double *x,double *v);
double rho_iso(double r);
double phi_iso(double r,double b);

#endif /* ISODF4_H_ */
