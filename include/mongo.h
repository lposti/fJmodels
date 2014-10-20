#include <math.h>
#include <stdio.h>
#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define SIGN(A,B) ((B)>0?(A):(-A))
#define MOD(A,B) ((A)-(B)*((int)(A)/(B)))

void plots(int,float*,float*,float,float,float,float,const char*,int,const char*,int,float,int);
void plots(int,double*,double*,float,float,float,float,const char*,int,const char*,int,float,int);
void grend(int);
void compress(FILE *,float *,int);
void compress(FILE *,double *,int);
void compress(FILE *,double **,int,int);
void compress(FILE *,double ***,int,int,int);
void compres2(FILE *,float *,int,int);
void compres2(FILE *,double *,int,int);
void compres2(FILE *,double **,int,int);
void get(FILE *,float *,int);
void get(FILE *,double *,int);
void get(FILE *,double **,int,int);
void get(FILE *,double ***,int,int,int);
void get2(FILE *,float *,int,int);
void get2(FILE *,double *,int,int);
void get2(FILE *,double **,int,int);
void choose(int,float *,float *,float *,float *);
void ctype(const char*);
void setcolour(const char*);
void setltype(int);
void setlweight(int);
void relocate(float,float);
void relocate(float*);
void label(int,const char *);
void putlabel(int,int,const char *);
void setlim(float,float,float,float);
void window(int,int,int);
void setangle(float);
void point(float);
void points(float,int,float *,float *,int);
void points(double,int,double *,double *,int);
void setexpand(float);
void draw(float,float);
void draw(float*);
void connect(float *,float *,int);
void connect(double *,double *,int);
void mgobox(int,int,int,int);
void xlabel(int,const char *);
void ylabel(int,const char *);
void errorbar(int,float *,float *,float *,int);
void errorbar(int,double *,double *,double *,int);
void contour(float *,int,int,float *,int);
void contour(double **,int,int,double *,int);
void contour(float **,int,int,float *,int);
void scontour(float *,int,int,float *,int,float,float,float,int);
void histogram(float *,float *,int);
void histogram(double *,double *,int);
void plt3d(float *,int,int,float,float,float,float);
void plt3d(double **,int,int,double,double,double,double);
void label3d(const char *,int,const char *,int,const char *,int);
void arrows(float *,int,float);
//void arrows(float *,int);
void greyscal(float *,int,int,float,float);
void greyscal(double **,int,int,double,double);
void coloursc(float *,int,int,float,float);
void coloursc(double **,int,int,double,double);
void colourwh(float *,int,int,float,float,float);
void colourwh(float **,int,int,float,float,float);
void colourwh(double *,int,int,double,double,double);
void colourwR(float *,int,int,float,float,float);
void colourwR(float**,int,int,double,double,double);
void colourwR(double **,int,int,double,double,double);
void shade(int,float *,float *,int);
void fillin(void);
void rgbcolor(int,int,int);
float plhist(float,float *,float *,int);
float plloghist(float,float *,float *,int);
