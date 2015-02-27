/*
 * writeOut.h
 *
 *  Created on: Feb 27, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_WRITEOUT_H_
#define INCLUDE_WRITEOUT_H_

#include <string>
#include <sstream>
#include <iomanip>
#include "readParam.h"

void writeOut(const struct fJParams, const int);


template <typename T> std::string toString( const T& t, int eps=2){

	std::stringstream ss;
	ss << std::fixed;				// format: could be also scientific
	ss << std::setprecision(eps);
	ss << t;
	return ss.str();
}


#endif /* INCLUDE_WRITEOUT_H_ */
