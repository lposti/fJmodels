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


string composeName(const struct fJParams fJP, const int iter){

	string name;

	for (int i=0; i<5; i++)
		name+=fJP.modName[i];

	name+="_"; name+=toString(fJP.dphi_h_in);
	name+="_"; name+=toString(fJP.dz_h_in);
	name+="_"; name+=toString(fJP.dphi_g_in);
	name+="_"; name+=toString(fJP.dz_g_in);
	name+="_"; name+=toString(iter);
	name+=".out";

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

void writeOut(const struct fJParams fJP, const int iter){

	ofstream outF;
	string name="models/";
	name+=composeName(fJP,iter);

	outF.open(name.c_str(), ios::out);

	if (outF.is_open()){

		outF << NR << " " << NPOLY << " " << NGAUSS << endl;

		writeArr(outF,ar,NR);
		writeMat(outF,rhl,NR,NPOLY);   writeMat(outF,vrotl,NR,NPOLY);
		writeMat(outF,sigRl,NR,NPOLY); writeMat(outF,sigpl,NR,NPOLY);
		writeMat(outF,sigzl,NR,NPOLY); writeMat(outF,sigRzl,NR,NPOLY);
	}

	outF.close();
}



