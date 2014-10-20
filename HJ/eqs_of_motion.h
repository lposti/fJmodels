/*
 * eqs_of_motion.h
 *
 *  Created on: Mar 16, 2013
 *      Author: Fermani
 */

#ifndef EQS_OF_MOTION_H_
#define EQS_OF_MOTION_H_

#include <stdio.h>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "Structure.h"
#include "../leg_pot3.h"

inline int eqs_of_motion(double t,const double y[], double f[],void *parameters){
	eveq_p *par= static_cast<eveq_p*>(parameters);
	f[0] = y[2];
	f[1] = y[3];
	if(par->Lz!=0) f[2] = -par->p->dphidR(y[0],y[1])+pow(par->Lz,2.0)/pow(y[0],3.0);
	else f[2] = -par->p->dphidR(y[0],y[1]);
	f[3] = -par->p->dphidz(y[0],y[1]);
	return GSL_SUCCESS;
}


inline int eqs_of_motion_spherical(double t,const double y[], double f[],void *parameters){
	eveq_p *par= static_cast<eveq_p*>(parameters);
	double d2;
	f[0] = y[2];
	f[1] = y[3]/pow(y[0],2);
	//std::cout << dPhitheta(y[0],cos(y[1])) << std::endl;
	f[2] = -dPhir(y[0],y[1],&d2)+pow(par->L,2.0)/pow(y[0],3.0);
	f[3] = -dPhitheta(y[0],cos(y[1]))+pow(par->Lz/(y[0]*sin(y[1])),2)/tan(y[1]);
	return GSL_SUCCESS;
}

#endif /* EQS_OF_MOTION_H_ */
