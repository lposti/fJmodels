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

extern double ** __restrict rhl, ** __restrict vrotl, ** __restrict sigRl,
              ** __restrict sigpl, ** __restrict sigzl, ** __restrict sigRzl;


#endif /* INCLUDE_GRID_H_ */
