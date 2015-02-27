/*
 * models.cpp
 *
 *  Created on: Feb 15, 2015
 *      Author: morpheus
 */
#include <math.h>
#include "models.h"
#include "Utils.h"

double rhoHern(double R, double z){
	double m=sqrt(R*R+z*z/q2);
	return mass/(m/r0*pow(r0+m,3))/TPI;
}

double rhoIsoch(double R, double z){
	double m=sqrt(R*R+z*z/q2);
	double a=sqrt(r0*r0+m*m);
	return mass/FPI/(q*pow((a+r0)*a,2)*a)*r0*(r0+2*a);
}


