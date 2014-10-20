//-----------------------------------------------------------------------------+
//                                                                             |
// Inline.h                                                                    |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// written by Walter Dehnen, 1994-2001                                         |
// e-mail:   dehnen@mpia-hd.mpg.de                                             |
// address:  Max-Planck Institut für Astronomie, Königstuhl 69117 Heidelberg   |
//           Germany                                                           |
//                                                                             |
//-----------------------------------------------------------------------------+

#ifndef _Inline_def_
#define _Inline_def_ 1

#define SUNS 1

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>

using std::exit;
using std::cerr;
////////////////////////////////////////////////////////////////////////////////
#ifndef ebug
void Numerics_error(const char*) __attribute__ ((noreturn));
#else
void Numerics_error(const char*);
#endif

inline void Numerics_error(const char* msgs)
{
    cerr << " ERROR in Numerics: " << msgs << '\n';
#ifndef ebug
    exit(1);
#endif
}

////////////////////////////////////////////////////////////////////////////////
inline int Numerics_message(const char* msgs)
{
  cerr << " WARNING in Numerics: " << msgs << '\n';
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

#include "inline2.h"

template<class S> inline S sqr(const S x)              { return square(x); }
template<class S> inline S pow(const S x, const int i) { return pow_int(x,i); }

////////////////////////////////////////////////////////////////////////////////
inline bool is_integral(const float x)
{
  if(x<0) return is_integral(-x);
  return (floor(x)==x)? 1 : 0; 
}

inline bool is_integral(const double x)
{
  if(x<0) return is_integral(-x);
  return (floor(x)==x)? 1 : 0; 
}

////////////////////////////////////////////////////////////////////////////////
#endif
