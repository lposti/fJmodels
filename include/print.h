/*
 * print.h
 *
 *  Created on: Oct 22, 2014
 *      Author: morpheus
 */

#ifndef PRINT_H_
#define PRINT_H_


//#define PRINTDFH
//#define PRINTXVJ


#ifdef PRINTDFH
void openDFH();
void printDFH(double Jr, double Jphi, double Jz, double DF, double E);
void closeDFH();
#endif

#ifdef PRINTXVJ
void openXVJ();
void printXVJ(double *x, double *v, double Jr, double Jphi, double Jz);
void closeXVJ();
#endif


#endif /* PRINT_H_ */
