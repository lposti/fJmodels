#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))

double bessj0(double);
double bessj1(double);
double bessi0(double);
double bessi1(double);
double bessk0(double);
double bessk1(double);

double qsimp(double (*)(double),double,double);
double qsimpi(double (*)(double),double,double);
int topbottom(double*,int,double,int*,int*);
int topbottom(float*,int,float,int*,int*);
float linterp(float *,float *,int,float);
double linterp(double *,double *,int,double);
double dlinterp(double *,double *,int,double);
double linterp2d(double **,double *,int,double *,int,double,double);
float linterp2d(float **,float *,int,float *,int,float,float);
double linterp3d(double***,double*,int,double*,int,double*,int,double,double,double);
void linterp2(double *,double **,int,double,double *);
void spline(double *,double *,int,double,double,double *);
double splint(double *, double *, double *,int,double);
double splintp(double *,double *,double *,int,double);
void rk4(double *,double *,int,double,double, void(*)(double,double *,double *));
void rkqs(double *,double *,int,double *,double,double,
	  double *,double *, double *,void (*)(double, double*, double*));
void rkqs(double*,double *,double *,int,double *,double,double,
	  double *,double *, double *,void (*)(double*,double, double*, double*));
double zbrent(double (*)(double),double,double,double,double,double,int);
double zbrent(double*,double (*)(double*,double),double,double,double,double,double,int);
double dbrent(double,double,double,double (*)(double),double (*)(double),double,double*);
int zbrac(double (*)(double),double *,double *,double *,double *,double);
double gammln(double);
double erff(double);
double ran1(long *);
double gasdev(long *);
void sort(unsigned long, double*);
void indexx(unsigned long,double*,unsigned long*);
double probks(double);
double ksone(double*,unsigned long,double (*)(double),double*);

void mrqmin(double*,double*,double*,int,double*,int*,int,double**,
	    double**,double*,void(*)(double,double*,double*,double*,int),double*);
void frprmn(double*,int,double,int*,double*,double (*)(double*),void (*)(double*,double*));
void gaussj(double**, int, double**, int);
void gauher(double*,double*,int);
void gauleg(double,double,double*,double*,int);
void tridiag(double*,double*,int,double**);

void chebft(double,double,double*,int,double (*)(double));
double chebev(double,double,double*,int,double);
double polint(double*,double*,int,double,double*);
void polyfit(double,double,double*,int,double*,double*,int);

void four1(float*,unsigned long,int);
void four1(double*,unsigned long,int);
void fourn(float*,unsigned long*,int,int);
void fourn(double*,unsigned long*,int,int);

int amoeba(double**,double*,int,double,double (*)(double*,int));

double **matrix(long,long);
char **smatrix(long,long);
double *dmatrix(int);
double **dmatrix(int,int);
double ***dmatrix(int,int,int);
double ****dmatrix(int,int,int,int);
void delmatrix(double**,int);
void delmatrix(double***,int,int);
void delmatrix(double****,int,int,int);
float *fmatrix(int);
float **fmatrix(int,int);
float ***fmatrix(int,int,int);
float ****fmatrix(int,int,int,int);
void delmatrix(float**,int);
void delmatrix(float***,int,int);
void delmatrix(float****,int,int,int);
