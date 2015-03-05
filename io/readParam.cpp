/*
 * readParam.cpp
 *
 *  Created on: Feb 24, 2015
 *      Author: morpheus
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "readParam.h"

using namespace std;

string readStr(stringstream& ss){
	string str;
	ss >> str;
	return str;
}

double readVal(stringstream& ss){
	double v;
	ss >> v;
	return v;
}

struct fJParams readParam(){

	/* initialize ifstram for parameterFile */
	ifstream pF("param.txt");
	string line;
	struct fJParams fJP;
	fJP.modName  = "null";
	fJP.modName2 = "null";

	/* check is open */
	if (pF.is_open()){
		/*
		 * read line until end-of-file is reached
		 * and ignore lines starting with hash symbols,
		 * whitespaces, tabulations and empty lines.
		 */
		while (!getline(pF,line).eof())
			if (line[0]!='#'  and
				line[0]!=' '  and
				line[0]!='\0' and
				line[0]!='\t'){

					/* convert string to stream */
					stringstream lineStream;
					lineStream.str(line);

					/* now use stream with >> */
					string par;
					lineStream >> par;

					if (par=="model")       {fJP.modName   = readStr(lineStream); fJP.comp=1;}
					if (par=="h(J):dphi")   fJP.dphi_h_in = readVal(lineStream);
					if (par=="h(J):dz")	    fJP.dz_h_in   = readVal(lineStream);
					if (par=="g(J):dphi")   fJP.dphi_g_in = readVal(lineStream);
					if (par=="g(J):dz")     fJP.dz_g_in   = readVal(lineStream);
					if (par=="mass")	    fJP.mass	  = readVal(lineStream);
					if (par=="r0")     		fJP.r0		  = readVal(lineStream);
					if (par=="q")     		fJP.q		  = readVal(lineStream);

					if (par=="2:model")       {fJP.modName2  = readStr(lineStream); fJP.comp=2;}
					if (par=="2:h(J):dphi")   fJP.dphi_h_in2 = readVal(lineStream);
					if (par=="2:h(J):dz")     fJP.dz_h_in2   = readVal(lineStream);
					if (par=="2:g(J):dphi")   fJP.dphi_g_in2 = readVal(lineStream);
					if (par=="2:g(J):dz")     fJP.dz_g_in2   = readVal(lineStream);
					if (par=="2:mass")        fJP.mass_2	 = readVal(lineStream);
					if (par=="2:r0")          fJP.r0_2		 = readVal(lineStream);
					if (par=="2:q")           fJP.q_2		 = readVal(lineStream);
			}

	}
	return fJP;
}


void printParam(struct fJParams fJP){

	cout << "\n======================================================================" << endl;
	cout <<   "= f(J) Model computation" << endl;
	cout <<   "=\n" << endl;
	cout <<   "PARAMETERS of the model:\n" << endl;
	cout <<   " - Classical model:\t" << fJP.modName << endl;
	cout <<   " - h(J) :\t\tJr + " << fJP.dphi_h_in << " Jphi + " << fJP.dz_h_in << " Jz" << endl;
	cout <<   " - g(J) :\t\tJr + " << fJP.dphi_g_in << " Jphi + " << fJP.dz_g_in << " Jz" << endl;
	cout <<   " - Mass :\t\t" << fJP.mass << endl;
	cout <<   " - r0   :\t\t" << fJP.r0 << endl;
	cout <<   " - Flattening :\t\t" << fJP.q << endl;
	cout << endl;

	if (fJP.modName2 != "null"){
			cout <<   "PARAMETERS of the second component model:\n" << endl;
			cout <<   " - Classical model:\t" << fJP.modName2 << endl;
			cout <<   " - h(J) :\t\tJr + " << fJP.dphi_h_in2 << " Jphi + " << fJP.dz_h_in2 << " Jz" << endl;
			cout <<   " - g(J) :\t\tJr + " << fJP.dphi_g_in2 << " Jphi + " << fJP.dz_g_in2 << " Jz" << endl;
			cout <<   " - Mass :\t\t" << fJP.mass_2 << endl;
			cout <<   " - r0   :\t\t" << fJP.r0_2 << endl;
			cout <<   " - Flattening :\t\t" << fJP.q_2 << endl;
			cout << endl;
	}
}
