/*
 * readParam.h
 *
 *  Created on: Feb 24, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_READPARAM_H_
#define INCLUDE_READPARAM_H_

#include <cstring>

struct fJParams { std::string modName; double dphi_h_in, dz_h_in,dphi_g_in, dz_g_in;};

struct fJParams readParam();


#endif /* INCLUDE_READPARAM_H_ */
