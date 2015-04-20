/*
 * readParam.h
 *
 *  Created on: Feb 24, 2015
 *      Author: L. Posti
 */

#ifndef INCLUDE_READPARAM_H_
#define INCLUDE_READPARAM_H_

#include <string>

struct fJParams {

	/* number of iterations */
	int itermax;

	/* first component */
	std::string modName, outName;
	double dphi_h_in,dz_h_in,dphi_g_in,dz_g_in;
	double chi,mass,r0,q;
	double alpha,beta;

	/* second component */
	std::string modName2, outName2;
	bool rot_2;
	double dphi_h_in2,dz_h_in2,dphi_g_in2,dz_g_in2;
	double chi_2,mass_2,r0_2,q_2;
	double alpha_2,beta_2;
};

struct fJParams readParam(char const * fName = "param.txt");
void checkParameterFile(const struct fJParams *);
void printParam(const struct fJParams *);


#endif /* INCLUDE_READPARAM_H_ */
