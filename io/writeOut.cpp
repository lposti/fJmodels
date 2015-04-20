/*
 * writeOut.cpp
 *
 *  Created on: Feb 27, 2015
 *      Author: morpheus
 */

#include <iostream>
#include <fstream>
#include <string>
#include "Grid.h"
#include "writeOut.h"
#include "readParam.h"

using namespace std;


string composeName(const struct fJParams fJP, const int iter, const unsigned comp=1){

	string name;

	if (comp==1)
	{
		// trim model name to 4 chars
		name = fJP.modName;
		if (name.size()>5) name.erase(5,name.size());

		name+="_"; name+=toString(fJP.dphi_h_in);
		name+="_"; name+=toString(fJP.dz_h_in);
		name+="_"; name+=toString(fJP.dphi_g_in);
		name+="_"; name+=toString(fJP.dz_g_in);

	}
	else if (comp==2)
	{
		// trim model name to 4 chars
		name = fJP.modName2;
		if (name.size()>5) name.erase(5,name.size());

		name+="_"; name+=toString(fJP.dphi_h_in2);
		name+="_"; name+=toString(fJP.dz_h_in2);
		name+="_"; name+=toString(fJP.dphi_g_in2);
		name+="_"; name+=toString(fJP.dz_g_in2);
	}

	name+="_"; name+=toString(iter);
	name+=".out";

	if (comp==2) name+="2c";

	return name;
}

void writeArr( ofstream & outF, double *arr, int NX ){

	for (int i=0; i<NX; i++)
		outF << arr[i] << endl;
}

void writeMat( ofstream & outF, double **mat, int NX, int NY ){

	for (int i=0; i<NX; i++){
		for (int j=0; j<NY; j++)
			outF << mat[i][j] << " ";

		outF << endl;
	}
}

void writeOut(const struct fJParams& fJP, const int iter,
		const unsigned comp,double **rhlH,double ** sigRlH,double **sigplH,
		double ** sigzlH,double **sigRzlH,double **vrotlH,double **philH,
		double **PrH, double **Pr2H){

	ofstream outF;

	if (comp==1)
	{
		string name="models/";

		/* if the user's specified outName is not null use that, else compose name */
		if (fJP.outName != "null") name += fJP.outName+"_"+toString(iter)+".out";
		else name+=composeName(fJP,iter);

		outF.open(name.c_str(), ios::out);

		if (outF.is_open()){

			/* Print a meaningful header first */
			outF << NR << " " << NPOLY << " " << NGAUSS << endl;
			outF << fJP.dphi_h_in << " " << fJP.dz_h_in << " " << fJP.dphi_g_in << " " << fJP.dz_g_in << endl;
			outF << fJP.chi << " " << fJP.mass << " " << fJP.r0 << endl;
			outF << fJP.alpha << " " << fJP.beta << endl;

			writeArr(outF,ar,NR);
			writeMat(outF,rhlH,NR,NPOLY);   writeMat(outF,vrotlH,NR,NPOLY);
			writeMat(outF,sigRlH,NR,NPOLY); writeMat(outF,sigplH,NR,NPOLY);
			writeMat(outF,sigzlH,NR,NPOLY); writeMat(outF,sigRzlH,NR,NPOLY);
			writeMat(outF,philH,NR,NPOLY);  writeMat(outF,PrH,NR,NPOLY);
			writeMat(outF,Pr2H,NR,NPOLY);
		}

	}
	else if (comp==2)
	{
		string name="models/";

		/* if the user's specified outName is not null use that, else compose name */
		if (fJP.outName2 != "null") name += fJP.outName2+"_"+toString(iter)+".out";
		else name+=composeName(fJP,iter,comp);

		outF.open(name.c_str(), ios::out);

		if (outF.is_open()){

			/* Print a meaningful header first */
			outF << NR << " " << NPOLY << " " << NGAUSS << endl;
			outF << fJP.dphi_h_in2 << " " << fJP.dz_h_in2 << " " << fJP.dphi_g_in2 << " " << fJP.dz_g_in2 << endl;
			outF << fJP.chi_2 << " " << fJP.mass_2 << " " << fJP.r0_2 << endl;
			outF << fJP.alpha_2 << " " << fJP.beta_2 << endl;

			writeArr(outF,ar,NR);
			writeMat(outF,rhlH,NR,NPOLY);   writeMat(outF,vrotlH,NR,NPOLY);
			writeMat(outF,sigRlH,NR,NPOLY); writeMat(outF,sigplH,NR,NPOLY);
			writeMat(outF,sigzlH,NR,NPOLY); writeMat(outF,sigRzlH,NR,NPOLY);
			writeMat(outF,philH,NR,NPOLY);  writeMat(outF,PrH,NR,NPOLY);
			writeMat(outF,Pr2H,NR,NPOLY);
		}

	}

	outF.close();
}



