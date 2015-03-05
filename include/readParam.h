/*
 * readParam.h
 *
 *  Created on: Feb 24, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_READPARAM_H_
#define INCLUDE_READPARAM_H_

#include <string>

struct fJParams {
	std::string modName;
	double dphi_h_in,dz_h_in,dphi_g_in,dz_g_in;
	double mass,r0,q;
	unsigned comp;

	std::string modName2;
	double dphi_h_in2,dz_h_in2,dphi_g_in2,dz_g_in2;
	double mass_2,r0_2,q_2;
};

struct fJParams readParam();
void printParam(struct fJParams);


#endif /* INCLUDE_READPARAM_H_ */
