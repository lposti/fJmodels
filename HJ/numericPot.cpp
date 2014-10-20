/*
 * numericPot.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: morpheus
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../press.h"
#include "../mongo.h"
#include "../leg_pot3.h"
#include "../potLEG.h"
#include "numericPot.h"

numericPot::numericPot(int npar,double *coeff,char *fname) : potential_obj(-1,npar,coeff,fname) {
	//readPotential();
	getCurrentPot();
	std::cout << "done reading potential \n";
}

numericPot::~numericPot() {
	delete &nr_; delete &ngauss_; delete &npoly_;
	delete &ar_[0]; delete &rhl_[0][0]; delete &phil_[0][0];
	delete &Pr_[0][0]; delete &Pr2_[0][0];
}

/*
 * Shadowing of base class member functions
 * Phi, dPhi/dR, dPhi/dz
 */
double numericPot::phi(double R,double z){
	//if (z<ar_[1]) z=ar_[1];
	//if (R<ar_[1]) R=ar_[1];
	//if (z>ar_[nr_-1]) z=ar_[nr_-1];
	//if (R>ar_[nr_-1]) R=ar_[nr_-1];
	//printf("== phi:%f\n",evpot(R*R,z));
	return evpot(R*R,z);
}

double numericPot::dphidR(double R, double z) {//derivative of Phi with respect to cylindrical R
	//if (z<ar_[1]) z=ar_[1];
	//if (R<ar_[1]) R=ar_[1];
	//return evfor_(R*R,z);
	//printf("== dphidR:%f\n",dPhiRnum(R,z));
	return dPhiRnum(R,z);
}

double numericPot::dphidz(double R, double z) {//derivative of Phi with respect to cylindrical z
	//if (z<ar_[1]) z=ar_[1];
	//if (R<ar_[1]) R=ar_[1];
	//return evfor_(R*R,z);
	//printf("== dphidz:%f\n",dPhiz(R,z));
	return dPhiz(R,z);
}

double numericPot::PhiR(double R){ return this->phi(R,0); }

double numericPot::dPhiR(double R){ return this->dphidR(R,0); }

double numericPot::dPhi(double *x2,double *f2){
	f2[0]=this->dphidR(x2[0],x2[1]);
	f2[1]=this->dphidz(x2[0],x2[1]);
	return this->phi(x2[0],x2[1]);
}

void numericPot::getCurrentPot(){

	nr_ = nr; npoly_ = npoly; ngauss_=ngauss;

	/* initialize potential arrays */
	ar_ = new double[nr_];
	Pr_=dmatrix(nr_,npoly_); Pr2_=dmatrix(nr_,npoly_);
	phil_=dmatrix(nr_,npoly_); rhl_=dmatrix(nr_,npoly_);

	for (int i=0; i<nr_; i++){
		ar_[i] = ar[i];
		for (int j=0; j<npoly_; j++){
			phil_[i][j] = phil[i][j];
			Pr_[i][j]   = Pr[i][j];
			Pr2_[i][j]  = Pr2[i][j];
			rhl_[i][j]  = rhl[i][j];
		}
	}
}

/*
 * numericPot::readPotential(), used to read the Legendre phl,Pr,Pr2 from the input file
 */
void numericPot::readPotential(){
	FILE *ifile;
	char fname[200];

	sprintf(fname,"/home/morpheus/Dropbox/Oxford2014/code_paper_thesis_06-01/qt-relax/fJmodels/models/Ex10_isoth_tanb_3.out");

	ifile=fopen(fname,"r");
	if(ifile==NULL){
		printf("cannot open %s\n",fname);
		exit(1);
	}
	fscanf(ifile,"%d %d %d",&nr_,&npoly_,&ngauss_);

	/* allocating and getting ar[nr], rhl[nr,npoly] from file */
	ar_ = new double[nr_]; rhl_=dmatrix(nr_,npoly_);
	get(ifile,ar_,nr_); get2(ifile,rhl_,nr_,npoly_);

	/* initialize potential arrays */
	Pr_=dmatrix(nr_,npoly_); Pr2_=dmatrix(nr_,npoly_);
	phil_=dmatrix(nr_,npoly_);

	this->potleg_();
	fclose(ifile);
}

double numericPot::evfor_(double R2,double z){//returns Phi from phil
	double pol[npoly_],dphip[npoly_];
	double c=z/sqrt(R2+z*z), r=sqrt(R2+z*z);
	intfor_(r,dphip); this->evenlegend_(pol,c);
	double dphi=dphip[0];
	for(int np=1; np<npoly_; np++) dphi+=dphip[np]*pol[np];
	return dphi;
}
/*
 * evaluatest the interpolated potential
 */
double numericPot::evpot_(double R2,double z){//returns Phi from phil
	double pol[npoly_],phip[npoly_];
	double c=z/sqrt(R2+z*z), r=sqrt(R2+z*z);
	intpo2_(r,phip); this->evenlegend_(pol,c);
	double phi=phip[0];
	for(int np=1; np<npoly_; np++) phi+=phip[np]*pol[np];
	return phi;
}

/*
 * numericPot::intpo2_, interpolates phil to r
 */
void numericPot::intpo2_(double r,double *phip){
	if(r>=ar_[nr_-1]){
		for(int k=0; k<npoly_; k++)
			phip[k]=phil_[nr_-1][k]*pow(ar_[nr_-1]/r,2*k+1);
	} else {
		int top,bot;
		topbottom(ar_,nr_,r,&bot,&top);
		double db=r-ar_[bot], f1=db/(ar_[top]-ar_[bot]);
		for(int k=0; k<npoly_; k++){// linear interpolation
			phip[k]=f1*phil_[top][k]+(1-f1)*phil_[bot][k];
		}
		if(top<10){// add quadratic contributions
			int thr;
			if(f1<0.5 && bot>0) thr=bot-1; else thr=top+1;
			double dt=r-ar_[top];
			double f2=dt*db/((ar_[thr]-ar_[top])*(ar_[thr]-ar_[bot])),
			f3=(ar_[thr]-ar_[bot])/(ar_[top]-ar_[bot]);
			for(int k=0; k<npoly_; k++){
				phip[k]+=f2*(phil_[thr][k]-phil_[bot][k]-f3*(phil_[top][k]-phil_[bot][k]));
			}
		}
	}
}

void numericPot::intfor_(double r,double *dphip){
	if(r>=ar_[nr_-1]){
		for(int k=0; k<npoly_; k++)
			dphip[k]=-(2*k+1)*phil_[nr_-1][k]*pow(ar_[nr_-1]/r,2*k+1)/r;
	} else {
		int top,bot;
		topbottom(ar_,nr_,r,&bot,&top);
		double db=r-ar_[bot], f1=db/(ar_[top]-ar_[bot]);
		for(int k=0; k<npoly_; k++){// linear interpolation
			dphip[k]=f1*Pr_[top][k]+(1-f1)*Pr_[bot][k];
		}
		if(top<10){// add quadratic contributions
			int thr;
			if(f1<0.5 && bot>0) thr=bot-1; else thr=top+1;
			double dt=r-ar_[top];
			double f2=dt*db/((ar_[thr]-ar_[top])*(ar_[thr]-ar_[bot])),
			f3=(ar_[thr]-ar_[bot])/(ar_[top]-ar_[bot]);
			for(int k=0; k<npoly_; k++){
				dphip[k]+=f2*(Pr_[thr][k]-Pr_[bot][k]-f3*(Pr_[top][k]-Pr_[bot][k]));
			}
		}
	}
}

/*
 * forces in spherical coordinates
 */
double numericPot::dPhir_(double r,double theta){
	double pol[npoly_],dphip[npoly_];
	double grad1=0;
	evenlegend_(pol,cos(theta));
	intfor_(r,dphip);
	for(int np=0; np<npoly_; np++){
		grad1+=dphip[np]*pol[np];
	}
	return grad1;
}
double numericPot::dPhitheta_(double r, double c) {
    //returns the first derivative of the whole potential wrt spherical theta
	double dphitheta = 0;
	double dpol[npoly_], phinterp[npoly_];
	intpo2_(r, phinterp);
	for(int k=0; k<npoly_; k++) dpol[k] = dlegend_(c, 2*k); //fills the array dpol with derivatives of the even legendre polys
	for(int j=0; j<npoly_; j++) dphitheta += phinterp[j]*dpol[j]; //this is the derivative wrt cos(theta), need to then multiply by -sin(theta)
	dphitheta *= -sqrt(1 - c*c); //im taking the positive root here - any issues with this?
	return dphitheta;
}

/*
 * numericPot::potleg_, computes the multipole expansion for the potential from that of rho
 */
int numericPot::potleg_(){

    double G=1,TPI=6.28318530718;
#define ISOTHERMAL

    double I_int[nr_][npoly_],I_ext[nr_][npoly_];
    int nc=nr_/12;// change integration technique at this grid pt
	for(int i=0;i<npoly_;i++){
		I_int[0][i]=0; I_ext[0][i]=0;//integrals from zero to zero
	}

#if defined(HERNQUIST) || defined(NFW)
	/*
	 * 	Corrections for singular models (e.g. Hernquist)
	 */
	for(int np=0; np<npoly_; np++){
		int l=2*np;
		double r=ar_[0],a,p1=0;
		/*
		 *  replacing these integrals with their singular r^-1 versions
		 */
		p1=pow(ar_[0],2*l+2); a=rhl_[0][np]/pow(ar_[0],l-1);
		for(int n=1; n<nc; n++){//leading terms in all integrals
			r=ar_[n]; double p2=pow(r,2*l+2);
			double b=rhl_[n][np]/pow(r,l-1), arg=.5*(a+b);
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+2);
			I_ext[n][np]=I_ext[n-1][np]+arg*(r-ar_[n-1]);
			a=b; p1=p2;
		}
		a*=pow(r,l); p1=pow(r,l+3);
		for(int n=nc; n<nr_; n++){
			r=ar_[n]; double p2=pow(r,l+3);
			double b=rhl_[n][np], arg=.5*(a+b);
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+3);
			if(l!=2) I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,2-l)-pow(ar_[n-1],2-l))/(double)(2-l);
			else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar_[n-1]));
			a=b; p1=p2;
		}
	}
#elif defined ISOTHERMAL

	for(int np=0; np<npoly_; np++){
			int l=2*np;
			double r,a,p1=0;
			if(np==0){//corrections for quadratic term in a at a<a[1]
				a=rhl_[0][np];
				I_int[1][np]=-.2*(a-rhl_[1][0])*pow(ar_[1],3);
				I_ext[1][np]=-.25*(a-rhl_[1][0])*pow(ar_[1],2);
			} else {//assume const density inside ar[1] and rhl_l~r^l
				a=rhl_[1][np]/pow(ar_[1],l);
			}
			for(int n=1; n<nc; n++){//leading terms in all integrals
				r=ar_[n]; double p2=pow(r,2*l+3);
				double b=rhl_[n][np]/pow(r,l), arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+3);
				I_ext[n][np]=I_ext[n-1][np]+.5*arg*(r*r-ar_[n-1]*ar_[n-1]);
				a=b; p1=p2;
			}
			// replaced for r^-2 at r --> + infty
			a=rhl_[nc-1][np]*pow(ar_[nc-1],2); p1=pow(ar_[nc-1],l+1);
			for(int n=nc; n<nr_; n++){
				r=ar_[n]; double p2=pow(r,l+1);
				double b=rhl_[n][np]*pow(r,2), arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+1);
				if(l!=0)I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,-l)-pow(ar_[n-1],-l))/(double)(-l);
				else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar_[n-1]));
				a=b; p1=p2;
			}
	}

#elif defined ISOCHRONE
	for(int np=0; np<npoly_; np++){
			int l=2*np;
			double r,a,p1=0;
			if(np==0){//corrections for quadratic term in a at a<a[1]
				a=rhl_[0][np];
				I_int[1][np]=-.2*(a-rhl_[1][0])*pow(ar_[1],3);
				I_ext[1][np]=-.25*(a-rhl_[1][0])*pow(ar_[1],2);
			} else {//assume const density inside ar[1] and rhl_l~r^l
				a=rhl_[1][np]/pow(ar_[1],l);
			}
			for(int n=1; n<nc; n++){//leading terms in all integrals
				r=ar_[n]; double p2=pow(r,2*l+3);
				double b=rhl_[n][np]/pow(r,l), arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+3);
				I_ext[n][np]=I_ext[n-1][np]+.5*arg*(r*r-ar_[n-1]*ar_[n-1]);
				a=b; p1=p2;
			}
			a*=pow(r,l); p1=pow(r,l+3);
			for(int n=nc; n<nr_; n++){
				r=ar_[n]; double p2=pow(r,l+3);
				double b=rhl_[n][np], arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+3);
				if(l!=2) I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,2-l)-pow(ar_[n-1],2-l))/(double)(2-l);
				else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar_[n-1]));
				a=b; p1=p2;
			}
		}
#endif
	for(int np=0; np<npoly_; np++){
		int l=2*np;
		if(l==0){
			//phil_[0][np]=-TPI*G*(I_ext[nr-1][np]-I_ext[0][np]);
			phil_[0][np]=-TPI*G*(I_int[0][np]/pow(ar_[0],l+1)+
						(I_ext[nr_-1][np]-I_ext[0][np])*pow(ar_[0],l));
			Pr2_[0][np]=TPI*G*rhl_[0][np];
		}else{
			phil_[0][np]=0;
			if(l==2) Pr2_[0][np]=-TPI*G*((l-1)*l*I_ext[nr_-1][np]-(2*l+1)*rhl_[0][np]);
			else Pr2_[0][np]=TPI*G*(2*l+1)*rhl_[0][np];
		}
		Pr_[0][np]=0;
		for(int n=1; n<nr_; n++){
			phil_[n][np]=-TPI*G*(I_int[n][np]/pow(ar_[n],l+1)+
					    (I_ext[nr_-1][np]-I_ext[n][np])*pow(ar_[n],l));
			Pr_[n][np]=-TPI*G*(-(l+1)*I_int[n][np]/pow(ar_[n],l+2)+
					  l*(I_ext[nr_-1][np]-I_ext[n][np])*pow(ar_[n],l-1));
			Pr2_[n][np]=-TPI*G*((l+2)*(l+1)*I_int[n][np]/pow(ar_[n],l+3)+
					   (l-1)*l*(I_ext[nr_-1][np]-I_ext[n][np])*pow(ar_[n],l-2)
					   -(2*l+1)*rhl_[n][np]);
		}
	}
    return 1.;
}


/*
 * Legendre polynomials
 */
void numericPot::legend_(double *allpol, double c, int np) {
    // evaluates the legendre polys up to l = npoly at c ------------
	np++;
	allpol[0] = 1; if(np<2) return;
	allpol[1] = c;
	for(int i = 2; i < np; i++)
		allpol[i] = ((2*i-1)*c*allpol[i-1] - (i-1)*allpol[i-2]) /(double)i;
}
double numericPot::dlegend_(double c, int n) {
    // evaluates the derivative of the legendre polynomial of order n
	if(n == 0) return 0;
	double allpol[n+1];
	legend_(allpol, c, n);
	return (n*allpol[n-1] - n*c*allpol[n]) / (1 - c*c);
}
