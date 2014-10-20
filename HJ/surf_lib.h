/*
 * surf_lib.h
 *
 *  Created on: Mar 12, 2013
 *      Author: Fermani
 */

#ifndef SURF_LIB_H_
#define SURF_LIB_H_

#include <stdio.h>
#include <string>
#include "Structure.h"

using std::string;

class surf_lib{
public:
	surf_lib(eveq_p par,string file,int length);
	~surf_lib();
	void build();
	void build(double minE, double maxE);
	void save(double **surf_LIB);
	int read(double **surf_LIB);

private:
	int length_;
	double **list_;
	eveq_p a_;
	string output_;
};


#endif /* SURF_LIB_H_ */
