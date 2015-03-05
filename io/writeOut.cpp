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

	if (comp==1){

		// trim model name to 4 chars
		name = fJP.modName;
		if (name.size()>5) name.erase(5,name.size());

		name+="_"; name+=toString(fJP.dphi_h_in);
		name+="_"; name+=toString(fJP.dz_h_in);
		name+="_"; name+=toString(fJP.dphi_g_in);
		name+="_"; name+=toString(fJP.dz_g_in);
		name+="_"; name+=toString(iter);
		name+=".out";

	} else if (comp==2) {

		// trim model name to 4 chars
		name = fJP.modName2;
		if (name.size()>5) name.erase(5,name.size());

		name+="_"; name+=toString(fJP.dphi_h_in2);
		name+="_"; name+=toString(fJP.dz_h_in2);
		name+="_"; name+=toString(fJP.dphi_g_in2);
		name+="_"; name+=toString(fJP.dz_g_in2);
		name+="_"; name+=toString(iter);
		name+=".out";
	}

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

void writeOut(const struct fJParams fJP, const int iter,
		const unsigned comp,double **rhlH,double ** sigRlH,double **sigplH,
		double ** sigzlH,double **sigRzlH,double **vrotlH,double **philH){

	ofstream outF;

	if (comp==1){
		string name="models/";
		name+=composeName(fJP,iter);

		outF.open(name.c_str(), ios::out);

		if (outF.is_open()){

			outF << NR << " " << NPOLY << " " << NGAUSS << endl;

			writeArr(outF,ar,NR);
			writeMat(outF,rhlH,NR,NPOLY);   writeMat(outF,vrotlH,NR,NPOLY);
			writeMat(outF,sigRlH,NR,NPOLY); writeMat(outF,sigplH,NR,NPOLY);
			writeMat(outF,sigzlH,NR,NPOLY); writeMat(outF,sigRzlH,NR,NPOLY);
			writeMat(outF,philH,NR,NPOLY);
		}

	} else {

		string name="models/";
		name +=composeName(fJP,iter,comp);

		outF.open(name.c_str(), ios::out);

		if (outF.is_open()){

			outF << NR << " " << NPOLY << " " << NGAUSS << endl;

			writeArr(outF,ar,NR);
			writeMat(outF,rhlH,NR,NPOLY);   writeMat(outF,vrotlH,NR,NPOLY);
			writeMat(outF,sigRlH,NR,NPOLY); writeMat(outF,sigplH,NR,NPOLY);
			writeMat(outF,sigzlH,NR,NPOLY); writeMat(outF,sigRzlH,NR,NPOLY);
			writeMat(outF,philH,NR,NPOLY);
		}

	}

	outF.close();
}



