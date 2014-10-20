//Potential & forces of flattened isochrone model. Must call init(b,q)
//first
#include <stdio.h>
#include <math.h>
#include "press.h"

#define npt 2000
#define SMALL .00001

double mx[npt],psix[npt],psix3[npt];
static double R,z,b,b2,q,q2,e,pi=3.1415926535897932,Mofpi=1./(4.*pi);

double rapzd(double (*func)(double), double a, double b,double *s,int *it, int n){
	double x,sum,del;
	int j;
	if (n == 0) {
		*s=0.5*(b-a)*((*func)(a)+(*func)(b));
		(*it)=1;
	} else {
		del=(b-a)/(*it);
		x=a+0.5*del; sum=0;
		for (j=0;j<(*it);j++,x+=del) sum += (*func)(x);
		*s=0.5*(*s+(b-a)*sum/(*it));
		(*it)*=2;
	}
	return *s;
}
#define EPS 1.0e-10
#define JMAX 20
double simp(double (*func)(double), double a, double b){
	int j,it;
	double s,st,ost,os;
	ost = os = -1.0e30;
	for (j=0;j<JMAX;j++) {
		s=(4.0*rapzd(func,a,b,&st,&it,j)-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		if (s == 0.0 && os == 0.0 && j > 6) return s;
		os=s;
		ost=st;
	}
	printf("Too many steps in routine qsimp\n");
	return s;
}
double simp(double (*func)(double), double a, double b,int *flag){
	int j,it;
	double s,st,ost,os;
	ost = os = -1.0e30; *flag=0;
	for (j=0;j<JMAX;j++) {
		s=(4.0*rapzd(func,a,b,&st,&it,j)-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		if (s == 0.0 && os == 0.0 && j > 6) return s;
		os=s;
		ost=st;
	}
	printf("Too many steps in routine qsimp\n");
	*flag=1; return s;
}
#undef EPS
#undef JMAX

double isoden(double m){
	double a=sqrt(b*b+m*m);
	return Mofpi/(q*pow((a+b)*a,2)*a)*b*(b+2*a);
}
double d_isoden_dm2(double m){
	double a=sqrt(b*b+m*m);
	//return -Mofpi/(q*pow((b+a)*a,3)*a)*b*(2*a+b)*(4*a+3*b)*.5/a;
	return -Mofpi/q*((2/pow((a+b)*a,3)+3/pow((a+b)*a*a,2))*b*(b+2*a)-2*b/(pow((a+b)*a,2)*a))*.5/a;
}
double isoden(double R,double z){
	double m=sqrt(R*R+z*z/q2);
	return isoden(m);
}
double misoden(double m){
	return m*isoden(m);
}
double chandra(double x){
	double psi,tau=1./(x*x)-b2;
	double m=sqrt(R*R/(tau/b2+1.)+z*z/(tau/b2+q2));
	if(m<mx[0]) psi=pow(m/mx[0],2)*psix[0];
	else if(m<mx[npt-1]) psi=splint(mx,psix,psix3,npt,m);
	else psi=psix[npt-1];
	return 2*psi*sqrt((tau+b2)/(tau+q2*b2));
}
double chandraR(double x){
	double tau=1./(x*x)-b2;
	double m=sqrt(R*R/(tau/b2+1.)+z*z/(tau/b2+q2));
	return -isoden(m)*4.*b2*R/(tau+b2)*sqrt((tau+b2)/(tau+q2*b2));
}
double chandraz(double x){
	double tau=1./(x*x)-b2;
	double m=sqrt(R*R/(tau/b2+1.)+z*z/(tau/b2+q2));
	return -isoden(m)*4.*b2*z/(tau+q2*b2)*sqrt((tau+b2)/(tau+q2*b2));
}
double chandraRR(double x){
	double tau=1./(x*x)-b2;
	double m=sqrt(R*R/(tau/b2+1.)+z*z/(tau/b2+q2));
	return -4.*b2/(tau+b2)*sqrt((tau+b2)/(tau+q2*b2))
			*(isoden(m)+d_isoden_dm2(m)*2*R*R/(tau/b2+1));
}
double chandrazz(double x){
	double tau=1./(x*x)-b2;
	double m=sqrt(R*R/(tau/b2+1.)+z*z/(tau/b2+q2));
	return -isoden(m)*4.*b2/(tau+q2*b2)*sqrt((tau+b2)/(tau+q2*b2));
}
void isopot_init(double bin,double qin){
	double small=0.01;
	q=qin; b=bin; q2=q*q; b2=b*b; e=sqrt(1-q2);
	double dm=50*b/float(npt-1);
	mx[0]=small; psix[0]=2*simp(&misoden,0,mx[0]);
	for(int i=1;i<npt;i++){
		mx[i]=mx[i-1]+dm;
		psix[i]=psix[i-1]+2*simp(&misoden,mx[i-1],mx[i]);
	}
	spline(mx,psix,npt,misoden(small),0,psix3);
}
double asinc(double e){
	if(e>1e-8) return asin(e)/e;
	else return 1+e*e/3;
}
double isopot(double Rin,double zin){
	R=Rin; z=zin;
	double msq=R*R+z*z/q2;
	if(msq==0) return -2.*pi*q*psix[npt-1]*asinc(e);
	double r2=R*R+z*z, m2=MIN(SMALL*msq,pow(mx[0],2));//pow(SMALL*mx[0],2);
	double bq=r2/m2-1.-q2, cq=q2-(R*R*q2+z*z)/m2;
	double taumax=b2*.5*(bq+sqrt(bq*bq-4.*cq));
	double xmin=1./sqrt(taumax+b2), xmax=1./b;
	int flag;
	double phi=simp(&chandra,xmin,xmax,&flag);
	if(flag==1)printf("%g %f %f %f %f\n",R,z,xmin,xmax,phi);
	return -2.*pi*q*(psix[npt-1]*asinc(e)-.5*b*phi);
}
void isoforce(double Rin,double zin,double *FR,double *Fz){
	R=Rin; z=zin;
	double msq=R*R+z*z/q2;
	if(msq==0){
		(*FR)=0; (*Fz)=0; return;
	}
	double r2=R*R+z*z, m2=MIN(SMALL*msq,pow(mx[0],2));//pow(SMALL*mx[0],2);
	double bq=r2/m2-1-q2, cq=q2-(R*R*q2+z*z)/m2;
	double taumax=b2*.5*(bq+sqrt(bq*bq-4*cq));
	double xmin=1/sqrt(taumax+b2), xmax=1/b;
	if(R==0) (*FR)=0;
	else (*FR)=pi*q*b*simp(&chandraR,xmin,xmax);
	if(z==0) (*Fz)=0;
	else (*Fz)=pi*q*b*simp(&chandraz,xmin,xmax);
}
void iso_d2Phi(double Rin,double *d2PhiR,double *d2Phiz){
	R=Rin; z=0;
	if(R==0) R=SMALL;
	double msq=R*R;
	double r2=R*R, m2=MIN(SMALL*msq,pow(mx[0],2));//pow(SMALL*mx[0],2);
	double bq=r2/m2-1-q2, cq=q2-(R*R*q2)/m2;
	double taumax=b2*.5*(bq+sqrt(bq*bq-4*cq));
	double xmin=1/sqrt(taumax+b2), xmax=1/b;
	(*d2PhiR)= -pi*q*b*simp(&chandraRR,xmin,xmax);
	(*d2Phiz)=-pi*q*b*simp(&chandrazz,xmin,xmax);
}
/*double d2Phiz(double Rin){
	R=Rin; z=0;
	double msq=R*R;
	if(msq==0) return 1;
	double r2=R*R, m2=MIN(SMALL*msq,pow(mx[0],2));//pow(SMALL*mx[0],2);
	double bq=r2/m2-1-q2, cq=q2-(R*R*q2)/m2;
	double taumax=b2*.5*(bq+sqrt(bq*bq-4*cq));
	double xmin=1/sqrt(taumax+b2), xmax=1/b;
	return -pi*q*b*simp(&chandrazz,xmin,xmax);
}*/
/*
int main(void){
	double R0=1,z0,FR,Fz,FR1,Fz1,FR2,Fz2,FR3,Fz3,FR4,Fz4;
	b=1; q=.7;
	isopot_init(b,q);
	while(R0>=0){
		printf("Enter R,z ");
		scanf("%lf %lf",&R0,&z0);
		if(R0<0) return 1;
		double dR=.001*R0, dz=.001*fabs(z0)+.00001;
		FR=-.5*(isopot(R0+dR,z0)-isopot(R0-dR,z0))/dR;
		Fz=-.5*(isopot(R0,z0+dz)-isopot(R0,z0-dz))/dz;
		isoforce(R0,z0,&FR1,&Fz1);
		isoforce(R0+dR,z0,&FR3,&Fz4); isoforce(R0-dR,z0,&FR4,&Fz4);
		FR3=-.5*(FR3-FR4)/dR;
		isoforce(R0,z0+dz,&FR4,&Fz3); isoforce(R0,z0-dz,&FR4,&Fz4);
		Fz3=-.5*(Fz3-Fz4)/dz;
		FR2=d2PhiR(R); Fz2=d2Phiz(R);
		double phi2=-1/(b+sqrt(b*b+R*R+z*z));
		printf("(%f %f) (%f %g) (%f %g) (%f %g) (%f %g)\n",
		       phi2,isopot(R0,z0),FR1,1-FR1/FR,Fz1,1-Fz1/Fz,
		       FR2,1-FR3/FR2,Fz2,1-Fz3/Fz2);
		getfreqs(R,&FR,&FR1,&FR2);
		printf("kappa, nu, omega: %f %f %f\n",FR,FR1,FR2);
	}
	return 1;
}
*/



/*
 *  Hernquist stuff
 */
#define TPI 6.283185307179586

double rho_Hern(double r){
	return 1./(r/b*pow(1.+r/b,3))/TPI;
}

double rho_Hern(double R, double z){
	double m=sqrt(R*R+z*z/q2);
	return rho_Hern(m);
}

double phi_Hern(double r){
	return -1./(1.+r/b);
}

double phi_iso(double r){
	return -1./(1.+sqrt(1.+r*r));
}


/*
 *  NFW stuff
 */
#define FPI 12.5663706143592

double rho_NFW(double r){
	return 1./(r/b*pow(1.+r/b,2))/FPI;
}
double rho_NFW(double R,double z){
	double m=sqrt(R*R+z*z/q2);
	return rho_NFW(m);
}

double phi_NFW(double r){
	return -log(1.+r/b)/(r/b);
}


/*
 *  Isothermal sphere stuff
 */
double sig=1.;
/*
 *  non-singular isothermal rho(r)
 *  normalized so that the Mass inside Rmax is
 *  M(r<Rmax)=1
 */
double rho_isoth(double r){
	//double mass=2*sig*sig*50.*b;
	return sig*sig/pow(1.+r,2)/TPI;// /mass;
}
double rho_isoth(double R, double z){
	double m=sqrt(R*R+z*z/q2);
	return rho_isoth(m);
}
/*
 *  logarithmic potential
 */
double phi_isoth(double r){
	return 2.*sig*sig*(log(r+b)-log((50.+1.)*b));
}
