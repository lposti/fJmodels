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
extern double ar[NR], rhl[NR][NPOLY], phil[NR][NPOLY], Pr[NR][NPOLY], Pr2[NR][NPOLY];


#endif /* INCLUDE_GRID_H_ */
