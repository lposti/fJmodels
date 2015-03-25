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
extern double ar[NR], Dgrid[NGRID], Egrid[NGRID];
extern double ** __restrict phil, ** __restrict Pr, ** __restrict Pr2,
			  ** __restrict rhl, ** __restrict vrotl, ** __restrict sigRl,
              ** __restrict sigpl, ** __restrict sigzl, ** __restrict sigRzl;

extern unsigned components;


#endif /* INCLUDE_GRID_H_ */
