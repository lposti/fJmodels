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
#include <stdlib.h>
#include "Grid.h"
#include "readParam.h"

#define UNSET -999.

unsigned components;
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

struct fJParams readParam(char const * fName){

	/* initialize ifstram for parameterFile */
	ifstream pF(fName);
	string line;
	struct fJParams fJP;

	/* default variables */
	fJP.itermax = 5;

	fJP.modName  = "null"; fJP.outName  = "null";
	fJP.dphi_h_in = UNSET; fJP.dz_h_in = UNSET; fJP.dphi_g_in = UNSET; fJP.dz_g_in = UNSET;
	fJP.chi = UNSET; fJP.mass = UNSET; fJP.r0 = UNSET; fJP.q = UNSET;

	fJP.modName2 = "null"; fJP.outName2  = "null";
	fJP.dphi_h_in2 = UNSET; fJP.dz_h_in2 = UNSET; fJP.dphi_g_in2 = UNSET; fJP.dz_g_in2 = UNSET;
	fJP.chi_2 = UNSET; fJP.mass_2 = UNSET; fJP.r0_2 = UNSET; fJP.q_2 = UNSET;


	/* check is open */
	if (pF.is_open())
	{
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

					if (par=="itermax")		fJP.itermax	  = readVal(lineStream);

					if (par=="model")       fJP.modName   = readStr(lineStream);
					if (par=="outfile")     fJP.outName   = readStr(lineStream);
					if (par=="chi")			fJP.chi		  = readVal(lineStream);
					if (par=="h(J):dphi")   fJP.dphi_h_in = readVal(lineStream);
					if (par=="h(J):dz")	    fJP.dz_h_in   = readVal(lineStream);
					if (par=="g(J):dphi")   fJP.dphi_g_in = readVal(lineStream);
					if (par=="g(J):dz")     fJP.dz_g_in   = readVal(lineStream);
					if (par=="mass")	    fJP.mass	  = readVal(lineStream);
					if (par=="r0")     		fJP.r0		  = readVal(lineStream);
					if (par=="q")     		fJP.q		  = readVal(lineStream);

					if (par=="2:model")       fJP.modName2   = readStr(lineStream);
					if (par=="2:outfile")     fJP.outName2   = readStr(lineStream);
					if (par=="2:chi")		  fJP.chi_2		 = readVal(lineStream);
					if (par=="2:h(J):dphi")   fJP.dphi_h_in2 = readVal(lineStream);
					if (par=="2:h(J):dz")     fJP.dz_h_in2   = readVal(lineStream);
					if (par=="2:g(J):dphi")   fJP.dphi_g_in2 = readVal(lineStream);
					if (par=="2:g(J):dz")     fJP.dz_g_in2   = readVal(lineStream);
					if (par=="2:mass")        fJP.mass_2	 = readVal(lineStream);
					if (par=="2:r0")          fJP.r0_2		 = readVal(lineStream);
					if (par=="2:q")           fJP.q_2		 = readVal(lineStream);
			}

	}
	else
	{
		cout << "ERROR in opening parameter file " << fName << "! Are you sure the file exists?" << endl;
		exit(5);
	}
	checkParameterFile(&fJP);
	return fJP;
}


void checkParameterFile(const struct fJParams * fJP){

	/* check for one component */
	if (fJP->modName != "null" and fJP->chi != UNSET and fJP->mass != UNSET and fJP->r0 != UNSET and fJP->q != UNSET
			and fJP->dphi_h_in != UNSET and fJP->dz_h_in != UNSET and fJP->dphi_g_in != UNSET and fJP->dz_g_in != UNSET)
	{
		/*  set the number of components  */
		components = 1;

		/* check for second component */
		if (fJP->modName2 != "null")
		{
			if (fJP->chi_2 != UNSET and fJP->mass_2 != UNSET and fJP->r0_2 != UNSET and fJP->q_2 != UNSET
				and fJP->dphi_h_in2 != UNSET and fJP->dz_h_in2 != UNSET and fJP->dphi_g_in2 != UNSET and fJP->dz_g_in2 != UNSET)

				components++;
			else
			{
				cout << "ERROR in parameterfile: second components is missing some parameter to initialize!" << endl;
				exit(2);
			}
		}

	}
	else
	{
		printf("ERROR: parameterfile must contain at least one component!\n");
		exit(1);
	}
}


void printParam(const struct fJParams * fJP){

	cout << "\n======================================================================" << endl;
	cout <<   "= f(J) Model computation" << endl;
	cout <<   "======================================================================\n" << endl;
	cout <<   "= Maximum number of iterations: " << fJP->itermax << endl;
	cout <<   "=\n" << endl;
	cout <<   "PARAMETERS of the model:\n" << endl;
	cout <<   " - Classical model:\t" << fJP->modName << endl;
	cout <<   " - h(J) :\t\tJr + " << fJP->dphi_h_in << " Jphi + " << fJP->dz_h_in << " Jz" << endl;
	cout <<   " - g(J) :\t\tJr + " << fJP->dphi_g_in << " Jphi + " << fJP->dz_g_in << " Jz" << endl;
	cout <<   " - Chi :\t\t" << fJP->chi << endl;
	cout <<   " - Mass :\t\t" << fJP->mass << endl;
	cout <<   " - r0   :\t\t" << fJP->r0 << endl;
	cout <<   " - Flattening :\t\t" << fJP->q << endl;
	cout << endl;
	if (fJP->outName != "null")
	{
		cout << " the model's output file will be " << fJP->outName << "_*.out" << endl;
		cout << endl;
	}

	if (components > 1){
			cout <<   "PARAMETERS of the second component model:\n" << endl;
			cout <<   " - Classical model:\t" << fJP->modName2 << endl;
			cout <<   " - h(J) :\t\tJr + " << fJP->dphi_h_in2 << " Jphi + " << fJP->dz_h_in2 << " Jz" << endl;
			cout <<   " - g(J) :\t\tJr + " << fJP->dphi_g_in2 << " Jphi + " << fJP->dz_g_in2 << " Jz" << endl;
			cout <<   " - Chi :\t\t" << fJP->chi_2 << endl;
			cout <<   " - Mass :\t\t" << fJP->mass_2 << endl;
			cout <<   " - r0   :\t\t" << fJP->r0_2 << endl;
			cout <<   " - Flattening :\t\t" << fJP->q_2 << endl;
			cout << endl;
			if (fJP->outName2 != "null")
			{
				cout << " the model's output file will be " << fJP->outName2 << "_*.out" << endl;
				cout << endl;
			}
	}
}
