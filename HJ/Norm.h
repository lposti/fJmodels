/*******************************************************************************
*                                                                              *
* Norm.h                                                                       *
*                                                                              *
* C++ code written by Walter Dehnen, 1996,                                     *
* Oxford University, Department of Physics, Theoretical Physics                *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail : dehnen@thphys.ox.ac.uk                                              *
*                                                                              *
*******************************************************************************/

#ifndef _Norm_def_
#define _Norm_def_ 1

inline int         norm(const int x) 		{ return x*x; }
inline long int    norm(const long int x) 	{ return x*x; }
inline float       norm(const float x) 		{ return x*x; }
inline double      norm(const double x) 	{ return x*x; }
inline long double norm(const long double x) 	{ return x*x; }

#endif
