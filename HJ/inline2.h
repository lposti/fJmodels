//-----------------------------------------------------------------------------+
//                                                                             |
// inline.h                                                                    |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// written by Walter Dehnen, 1994-2001                                         |
// e-mail:   dehnen@mpia-hd.mpg.de                                             |
// address:  Max-Planck Institut für Astronomie, Königstuhl 69117 Heidelberg   |
//           Germany                                                           |
//                                                                             |
//-----------------------------------------------------------------------------+

#ifndef WD_inline_h
#define WD_inline_h

#ifndef C_AND_F

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#else // C_AND_F
// C & FORTRAN versions: don't use C++ I/O, but plain C I/O.

#include <stdio>

#endif // C_AND_F

//using namespace std;

using std::exit;

using std::ofstream;
using std::ifstream;
using std::fstream;
using std::istream;
using std::ostream;
using std::ios;
using std::string;
using std::_Ios_Openmode;

//------------------------------------------------------------------------------
#define TM template<class S> inline

//------------------------------------------------------------------------------
TM S    abs     (const S x)              { return (x<0)? -x : x; }
TM int  sign    (const S x)              { return (x<0)? -1:((x>0)? 1:0 ); }
TM S    sign    (const S x, const S s)   { return ( s>0 )?  abs(x) : -abs(x); }
TM S    WDmin     (const S x, const S y)   { return (x<y)? x : y; }
TM S    WDmax     (const S x, const S y)   { return (x>y)? x : y; }
TM S    mod     (const S x, const S y)   { return x-y*int(x/y); }
TM S    square  (const S x)              { return x*x; }
TM S    cube    (const S x)              { return x*x*x; }
TM void WDswap    (S&a, S&b)               { register S t=a; a=b; b=t; }
TM S    pow_uint(const S x, const unsigned int n)
{
  if(n==0) return S(1);
  register S z=x, y=(n%2)? x : S(1);
  for(register unsigned int i=n>>1; i; i>>=1) { z*=z; if(i%2) y*=z; }
  return y;
}
TM S    pow_int (const S x, const int i)
{
  if(i>=0) return pow_uint(x,i);
#ifdef C_AND_F
  if(x==0) { printf("pow_int(): negative power of zero"); exit(1); }
#else // C_AND_F
  if(x==0) { cerr<<"pow_int(): negative power of zero"; exit(1); }
#endif // C_AND_F
  return pow_uint(S(1)/x,-i);
}

TM void update_max(S&x, const S y) { if(y>x) x=y; }
TM void update_min(S&x, const S y) { if(y<x) x=y; }

#undef TM

//------------------------------------------------------------------------------
template<class X>
inline char* negspace(const X x)
{
  return (x<X(0))? " " : "  ";
}
//------------------------------------------------------------------------------
#ifndef C_AND_F
//------------------------------------------------------------------------------
inline int open(ofstream& S, const char* file, _Ios_Openmode mode=ios::out)
{
    S.open(file,mode);
    if(! S.is_open() ) {
        cerr<<" cannot open file "<<file;
        return 0;
    }
    return 1;
}
//------------------------------------------------------------------------------
inline int open(ofstream& S, const string file, _Ios_Openmode mode=ios::out)
{
    return open(S,file.c_str(),mode);
}

//------------------------------------------------------------------------------
inline int open_to_append(ofstream& S, const char* file)
{
    S.open(file,ios::out | ios::app );
    if(! S.is_open() ) S.open(file,ios::out);
    if(! S.is_open() ) {
        cerr<<" cannot open file "<<file;
        return 0;
    }
    return 1;
}

//------------------------------------------------------------------------------
inline int open_to_append(ofstream& S, const string file)
{
    return open_to_append(S,file.c_str());
}

//------------------------------------------------------------------------------
inline int open(ifstream& S, const char* file, _Ios_Openmode mode=ios::in)
{
    S.open(file,mode);
    if(! S.is_open() ) {
        cerr<<" cannot open file "<<file;
        return 0;
    }
    return 1;
}

//------------------------------------------------------------------------------
inline int open(ifstream& S, const string file, _Ios_Openmode mode=ios::in)
{
   return open(S,file.c_str(),mode);
}

//------------------------------------------------------------------------------
inline int open(fstream& S, const char* file, _Ios_Openmode mode)
{
    S.open(file,mode);
    if(! S.is_open() ) {
        cerr<<" cannot open file "<<file;
        return 0;
    }
    return 1;
}

//------------------------------------------------------------------------------
inline int open(fstream& S, const string file, _Ios_Openmode mode)
{
    return open(S,file.c_str(),mode);
}

//------------------------------------------------------------------------------
inline void SwallowRestofLine(istream& from)
{
    char c;
    do from.get(c); while( from.good() && c !='\n');
}
//------------------------------------------------------------------------------
inline const char* stndrdth(const int i)
{
    register int ia= (i<0)? -i:i;
    switch( ia % 100 ) {
        case 11: 
        case 12: 
        case 13: return "th";
        default:
        switch( ia % 10 ) {
            case 1:  return "st";
            case 2:  return "nd"; 
            case 3:  return "rd";
            default: return "th";
        }   
    }
}
//------------------------------------------------------------------------------
#endif // C_AND_F
#endif // inline_h
