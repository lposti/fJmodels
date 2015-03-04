#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([fJmodels], [0.9], [lorenzo.posti@unibo.it])
AC_CONFIG_SRCDIR([fJmodels.cpp])

if test ["$CXX" == ""] || ["$CC" == ""] ; then
  CC=gcc
  CXX=g++
fi
if test ["$CXX" == g++] || ["$CC" == gcc] ; then
      CC=gcc
      CXX=g++
      FLAGS="-g3 -fopenmp"
      OPTFLAGS="-O3 -march=native -mtune=native -msse2"
      #OPTFLAGS="-O3 -ftree-vectorizer-verbose=1 -march=corei7-avx -mtune=corei7-avx -mavx"
fi
if test ["$CXX" == icpc] || ["$CC" == icc] ; then
      CC=icc
      CXX=icpc
      FLAGS="-g3 -parallel -openmp"
      OPTFLAGS="-O3 -qopt-report=3 -ip -xavx -axavx -fno-alias -fno-fnalias"
fi
AC_SUBST(FLAGS)
AC_SUBST(OPTFLAGS)

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)


# Extra places to look for third-party include files and libraries.
INCDIRS="$INCDIRS           \
         /usr/include       \
         /usr/local/include \
         /usr/include/x86_64-linux-gnu/c++/4.8 \
         /opt/local/include"

LIBDIRS="$LIBDIRS           \
         /usr/lib           \
         /usr/local/lib     \
         /local/lib         \
         /opt/local/lib"

for LIBDIR in $LIBDIRS ; do
  AC_CHECK_FILE([$LIBDIR], [LDFLAGS="$LDFLAGS -L$LIBDIR"], [continue])
done


# --with GSL library
AC_ARG_WITH(gslprefix,
  AC_HELP_STRING([--with-gslprefix=<dir>],[Specify custom GSL prefix]),
  [
    # user specified GSL prefix path
    if test "x$with_gslprefix" == xno
    then
        AC_MSG_NOTICE([The GSL library is required. Searching for the library in standard paths.], [1])
        GSLINCS=$INCDIRS
        GSLLIBS=$LIBDIRS
    else
        AC_MSG_NOTICE([Using GSL prefix: "$with_gslprefix"])
        GSLINCS="${with_gslprefix}include $with_gslprefix/include"
        GSLLIBS="${with_gslprefix}lib $with_gslprefix/lib"
    fi
    # including libs and headers in std. flags
    for LIBDIR in $GSLLIBS ; do
      AC_CHECK_FILE([$LIBDIR], [LDFLAGS="$LDFLAGS -L$LIBDIR"], [continue])
    done
    for INCDIR in $GSLINCS ; do
      AC_CHECK_FILE([$INCDIR/gsl],
                    [CPPFLAGS="$CPPFLAGS -I$INCDIR"; GSLINC=$INCDIR; break],
                    [continue])
    done
    #CPPFLAGS="$CPPFLAGS -I"$with_gslprefix"/include -L"$with_gslprefix"/lib"
  ],
  [
    for INCDIR in $INCDIRS ; do
      AC_CHECK_FILE([$INCDIR/gsl],
                    [CPPFLAGS="$CPPFLAGS -I$INCDIR"; GSLINC=$INCDIR; break],
                    [continue])
    done
    for LIBDIR in $GSLLIBS ; do
      AC_CHECK_FILE([$LIBDIR], [LDFLAGS="$LDFLAGS -L$LIBDIR"], [continue])
    done
  ]
)




# Checks for libraries.
AC_CHECK_LIB([m],[cos],[],
  [AC_MSG_ERROR([
     ========================================
     =      ERROR in libm lookup!           =
     ========================================],[1])])


# GSL standard library check
AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],
  [AC_MSG_ERROR([
     ==========================================
     =      ERROR in libgslcblas lookup!      =
     =                                        =
     = GSL library needed for compilation,    =
     = try running:                           =
     = ./configure --with-gslprefix=<GSLPATH> =
     ==========================================],[1])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],
  [AC_MSG_ERROR([
     ==========================================
     =      ERROR in libgslcblas lookup!      =
     =                                        =
     = GSL library needed for compilation,    =
     = try running:                           =
     = ./configure --with-gslprefix=<GSLPATH> =
     ==========================================],[1])])
# now $LIBS=-lgsl -lgslcblas -lm

# Checks for header files.
AC_CHECK_HEADERS(stdlib.h stdio.h math.h time.h string fstream iostream,[],
  [AC_MSG_ERROR([
    ==========================================
    =  ERROR in standard headers lookup!     =
    ==========================================
  ],[1])
  ])

# GSL headers check
AC_CHECK_HEADERS(gsl/gsl_blas.h gsl/gsl_cblas.h gsl/gsl_errno.h gsl/gsl_math.h gsl/gsl_integration.h gsl/gsl_roots.h,[],
  [AC_MSG_ERROR([
    ==========================================
    =    ERROR in GSL headers lookup!        =
    ==========================================
  ],[1])
  ])


# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt exp log log10 fabs clock],[],
  [AC_MSG_ERROR([
    ==========================================
    =  ERROR in standard functions lookup!   =
    ==========================================
  ],[1])
  ])


#AC_SUBST(LIBDIRS)
#AC_SUBST(INCDIRS)
#AC_SUBST(LIBDIR)
AC_SUBST(LDFLAGS)
AC_OUTPUT(Makefile)
