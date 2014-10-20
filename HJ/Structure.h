/*
 * Structure.h
 *
 *  Created on: Mar 12, 2013
 *      Author: Fermani
 */

#ifndef STRUCTURE_H_
#define STRUCTURE_H_

#include "potential_obj.h"

// eveq_p for SHO
/*struct eveq_p{
	double x_spring;
	double y_spring;
};*/

struct eveq_p{
	double energy;
	double L;
	double Lz;
	potential_obj *p;
};

struct data {
	size_t n;
	double *jr;
	double *jz;
	double *jphi;
	double *sigma;
};

const int maxT=1.0e4;
const int ID=1.0e6;
const int max_ite=20;
const double rsup=1.0e7;
const double rinf=1.0e-3;//-5
const double tollRF=1.0e-8;
const double toll=1.0e-3;
const double step=1.0;
const double unit=9.7777, unit2=unit*unit;

const double SMALL=1.0e-4;

const int lmax=12;
const int planeorquadr=0; //1 for quadratic fit, 0 for linear fti
const int nsurf=500;

#endif /* STRUCTURE_H_ */
