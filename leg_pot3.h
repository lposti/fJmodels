extern double *ar,**rhl,**vbarl,**sigRl,**sigpl,**sigzl,**sigRzl,**phil,**Pr,**Pr2;//defined in leg_pot2.cpp
extern double **phil_ini, **Pr_ini, **Pr2_ini; // initial potential used in ini_potLEG.cpp
extern int nr,ngauss,npoly;    //defined in leg_pot2.cpp
void setgrid(double,double);
void legend(double*,double,int);
double ev_dens(double,double);
double ev_dens(double,double,double*,double*,double*,double*,double*);
float ev_dens(float,float,float*,float*,float*,float*,float*);
double evpot(double,double);
void get_Ylm(double**,double**,int,int,int);
void potent(char*,double (*)(double,double),int,int);
void potent5(char*,double (*)(double,double,double*,double*,double*,double*,double*,double),int,int);
double phileg(double,double);
double dPhir(double,double,double*);
double dPhitheta(double,double);
double d2Phitheta(double,double);
void evenlegend(float *,float,int);
void evenlegend(double *,double,int);
void legend(double *,double,int);
double dlegend(double,int);
double d2legend(double,int);
