/*
 * Grid.h
 *
 *  Created on: Feb 13, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_GRID_H_
#define INCLUDE_GRID_H_

#define NR 		60
#define NGRID 	30
#define NGAUSS  6
#define NPOLY   4

/*
 *  Array definitions
 */
extern double ar[NR], phil[NR][NPOLY], Pr[NR][NPOLY], Pr2[NR][NPOLY];
extern double Dgrid[NGRID], Egrid[NGRID];

extern double **rhl, **vrotl, **sigRl, **sigpl, **sigzl, **sigRzl;


#endif /* INCLUDE_GRID_H_ */
