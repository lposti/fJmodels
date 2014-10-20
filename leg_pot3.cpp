#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "press.h"
#include "mongo.h"
#include "progressbar.hh"
#include "isodf4.h"
#include "potLEG.h"
#include "isopot.h"
#include "gsl/gsl_integration.h"

#define TINY 1e-6

//extern FILE *efile;
double *ar,**Pr,**Pr2,**rhl,**vbarl,**sigRl,**sigpl,**sigzl,**sigRzl,**phil;
int nr,ngauss,npoly;
const double G=1,TPI=2*acos(-1);
static int iter=0;

void setgrid(double a_scale,double rmax){
	double dpsi=asinh(rmax/a_scale)/(double)(nr-1);
	for(int n=0; n<nr; n++)	ar[n]=a_scale*sinh(n*dpsi);
	ar[0]+=.25*ar[1];													// adding tiny number to innermost grid point
}
void evenlegend(double *pol,double c,int npoly){
// evaluates even legendre polys up to l=2*(npoly-1) at c -------------
	double c2=c*c;
	pol[0]=1; if(npoly<2) return;
	pol[1]=1.5*c2-.5;
	for(int np=2; np<npoly; np++){
		int l=2*(np-1), l2=2*l;
		pol[np]=-pol[np-2]*l*(l-1)/(double)((l2+1)*(l2-1))+
			  pol[np-1]*(c2-(l2*l+l2-1)/(double)((l2-1)*(l2+3)));
		pol[np]*=(l2+1)*(l2+3)/(double)((l+1)*(l+2));
	}
}
void evenlegend(float *pol,float c,int npoly){
// evaluates even legendre polys up to l=2*(npoly-1) at c -------------
	float c2=c*c;
	pol[0]=1; if(npoly<2) return;
	pol[1]=1.5*c2-.5;
	for(int np=2; np<npoly; np++){
		int l=2*(np-1), l2=2*l;
		pol[np]=-pol[np-2]*l*(l-1)/(float)((l2+1)*(l2-1))+
			  pol[np-1]*(c2-(l2*l+l2-1)/(float)((l2-1)*(l2+3)));
		pol[np]*=(l2+1)*(l2+3)/(float)((l+1)*(l+2));
	}
}
void legend(double *allpol, double c, int npoly) {
    // evaluates the legendre polys up to l = npoly at c ------------
	npoly++;
	allpol[0] = 1; if(npoly<2) return;
	allpol[1] = c;
	for(int i = 2; i < npoly; i++)
		allpol[i] = ((2*i-1)*c*allpol[i-1] - (i-1)*allpol[i-2]) /(double)i;
}
double dlegend(double c, int n) {
    // evaluates the derivative of the legendre polynomial of order n
	if(n == 0) return 0;
	double allpol[n+1];
	legend(allpol, c, n);
	return (n*allpol[n-1] - n*c*allpol[n]) / (1 - c*c);
}
double d2legend(double c, int n) {
    //evaluates the second derivative of the legendre polynomial of order n
	double allpol[n+1];
	legend(allpol, c, n);
	return (2*c*dlegend(c, n) - n*(n+1)*allpol[n]) / (1. - c*c);
}
void int_dens(double r,double *densp,int npoly){
	if(r>ar[nr-1]){
		for(int i=0; i<npoly; i++){
			densp[i]=0;
		}
	}else{
		int top,bot;
		topbottom(ar,nr,r,&bot,&top);
		double f=(r-ar[bot])/(ar[top]-ar[bot]);
		for(int i=0; i<npoly; i++){
			densp[i]=f*rhl[top][i]+(1-f)*rhl[bot][i];
		}
	}
}
void int_dens(double r,double *densp,double *Vrotp,double *sigup,double *sigpp,double *sigvp,int npoly){
	if(r>ar[nr-1]){
		for(int i=0; i<npoly; i++){
			densp[i]=0; Vrotp[i]=0; sigup[i]=0; sigpp[i]=0; sigvp[i]=0;
		}
	}else{
		int top,bot;
		topbottom(ar,nr,r,&bot,&top);
		double f=(r-ar[bot])/(ar[top]-ar[bot]);
		for(int i=0; i<npoly; i++){
			densp[i]=f*rhl[top][i]+(1-f)*rhl[bot][i];
			sigup[i]=f*sigRl[top][i]+(1-f)*sigRl[bot][i];
			sigpp[i]=f*sigpl[top][i]+(1-f)*sigpl[bot][i];
			sigvp[i]=f*sigzl[top][i]+(1-f)*sigzl[bot][i];
			Vrotp[i]=f*vbarl[top][i]+(1-f)*vbarl[bot][i];
		}
	}
}
void int_dens(float r,float *densp,float *Vrotp,float *sigup,float *sigpp,float *sigvp,int npoly){
	if(r>ar[nr-1]){
		for(int i=0; i<npoly; i++){
			densp[i]=0; Vrotp[i]=0; sigup[i]=0; sigpp[i]=0; sigvp[i]=0;
		}
	}else{
		int top,bot;
		topbottom(ar,nr,r,&bot,&top);
		float f=(r-ar[bot])/(ar[top]-ar[bot]);
		for(int i=0; i<npoly; i++){
			densp[i]=f*rhl[top][i]+(1-f)*rhl[bot][i];
			sigup[i]=f*sigRl[top][i]+(1-f)*sigRl[bot][i];
			sigpp[i]=f*sigpl[top][i]+(1-f)*sigpl[bot][i];
			sigvp[i]=f*sigzl[top][i]+(1-f)*sigzl[bot][i];
			Vrotp[i]=f*vbarl[top][i]+(1-f)*vbarl[bot][i];
		}
	}
}
double ev_dens(double R,double z){//returns dens from rhl
	double pol[npoly],densp[npoly];
	double c=z/sqrt(R*R+z*z), r=sqrt(R*R+z*z);
	evenlegend(pol,c,npoly); int_dens(r,densp,npoly); 
	double dens=.5*densp[0];
	for(int np=1; np<npoly; np++){
		double fac=.5*(4*np+1);
		dens+=fac*densp[np]*pol[np];
	}
	return dens;
}
double ev_dens(double R,double z,double *Vrot,double *sigu,double *sigp,double *sigv,double *sigRz){//returns dens from rhl
	double pol[npoly],densp[npoly],Vrotp[npoly],sigup[npoly],sigpp[npoly],sigvp[npoly];
	double c=z/sqrt(R*R+z*z), r=sqrt(R*R+z*z);
	evenlegend(pol,c,npoly); int_dens(r,densp,Vrotp,sigup,sigpp,sigvp,npoly); 
	double dens=.5*densp[0];
	(*Vrot)=.5*Vrotp[0]; (*sigu)=.5*sigup[0]; (*sigp)=.5*sigpp[0]; (*sigv)=.5*sigvp[0];
	for(int np=1; np<npoly; np++){
		double fac=.5*(4*np+1);
		dens+=fac*densp[np]*pol[np]; (*sigu)+=fac*sigup[np]*pol[np];
		(*sigp)+=fac*sigpp[np]*pol[np]; (*sigv)+=fac*sigvp[np]*pol[np];
		(*Vrot)+=fac*Vrotp[np]*pol[np];
	}
	return dens;
}
float ev_dens(float R,float z,float *Vrot,float *sigu,float *sigp,float *sigv,float *sigRz){//returns dens from rhl
	float pol[npoly],densp[npoly],Vrotp[npoly],sigup[npoly],sigpp[npoly],sigvp[npoly];
	float c=z/sqrt(R*R+z*z), r=sqrt(R*R+z*z);
	evenlegend(pol,c,npoly); int_dens(r,densp,Vrotp,sigup,sigpp,sigvp,npoly);
	float dens=.5*densp[0];
	(*Vrot)=.5*Vrotp[0]; (*sigu)=.5*sigup[0]; (*sigp)=.5*sigpp[0]; (*sigv)=.5*sigvp[0];
	for(int np=1; np<npoly; np++){
		float fac=.5*(4*np+1);
		dens+=fac*densp[np]*pol[np]; (*sigu)+=fac*sigup[np]*pol[np];
		(*sigp)+=fac*sigpp[np]*pol[np]; (*sigv)+=fac*sigvp[np]*pol[np];
		(*Vrot)+=fac*Vrotp[np]*pol[np];
	}
	return dens;
}
void intpo2(double r,double *phip,int npoly){//interpolates phil to r
	if(r>=ar[nr-1]){
		for(int k=0; k<npoly; k++)
			phip[k]=phil[nr-1][k]*pow(ar[nr-1]/r,2*k+1);
		// useful for Debug, whether potential is positive
		//if (r==50.) printf(">> r=%f phip[0]=%f phip[1]=%f phil[%d][0]=%f phil[%d][1]=%f\n",r,phip[0],phip[1],nr-1,phil[nr-1][0],nr-1,phil[nr-1][1]);
	} else {
		/*
		 * 02/09/14: added check on minimum radius, so that I won't get seg. faults.
		 * 			 to be checked its consistency anyway
		 */
		int top,bot; r=MAX(r,ar[0]);
		topbottom(ar,nr,r,&bot,&top);
		double db=r-ar[bot], f1=db/(ar[top]-ar[bot]);
		for(int k=0; k<npoly; k++){// linear interpolation
			phip[k]=f1*phil[top][k]+(1-f1)*phil[bot][k];
		}
		if(top<10){// add quadratic contributions
			int thr;
			if(f1<0.5 && bot>0) thr=bot-1; else thr=top+1;
			double dt=r-ar[top];
			double f2=dt*db/((ar[thr]-ar[top])*(ar[thr]-ar[bot])),
			f3=(ar[thr]-ar[bot])/(ar[top]-ar[bot]);
			for(int k=0; k<npoly; k++){
				phip[k]+=f2*(phil[thr][k]-phil[bot][k]-f3*(phil[top][k]-phil[bot][k]));
			}
		}
		//if (r>1e-6) printf(">> r=%f phip[0]=%f phip[1]=%f phil[0]=%f phil[1]=%f\n",r,phip[0],phip[1],phil[bot][0],phil[bot][1]);
	}
}
void intpo2(double r,double *phip,double *dphip,double *d2phip,int npoly){//interpolates phil to r
	if(r>=ar[nr-1]){//r larger than end of grid
		for(int k=0; k<npoly; k++){
			phip[k]=phil[nr-1][k]*pow(ar[nr-1]/r,2*k+1);
			dphip[k]=-(2*k+1)*phil[nr-1][k]*pow(ar[nr-1]/r,2*k+1)/r;
			d2phip[k]=(2*k+2)*(2*k+1)*phil[nr-1][k]*pow(ar[nr-1]/r,2*k+1)/r/r;
		}
	} else {
		int top,bot;
		topbottom(ar,nr,r,&bot,&top);
		double db=r-ar[bot], f1=db/(ar[top]-ar[bot]);
		for(int k=0; k<npoly; k++){// linear interpolation
			phip[k]=f1*phil[top][k]+(1-f1)*phil[bot][k];
			dphip[k]=f1*Pr[top][k]+(1-f1)*Pr[bot][k];
			d2phip[k]=f1*Pr2[top][k]+(1-f1)*Pr2[bot][k];
		}
		if(top<10){// add quadratic contributions
			int thr;
			if(f1<0.5 && bot>0) thr=bot-1; else thr=top+1;
			double dt=r-ar[top];
			double f2=dt*db/((ar[thr]-ar[top])*(ar[thr]-ar[bot])),
			f3=(ar[thr]-ar[bot])/(ar[top]-ar[bot]);
			for(int k=0; k<npoly; k++){
				phip[k]+=f2*(phil[thr][k]-phil[bot][k]-f3*(phil[top][k]-phil[bot][k]));
				dphip[k]+=f2*(Pr[thr][k]-Pr[bot][k]-f3*(Pr[top][k]-Pr[bot][k]));
				d2phip[k]+=f2*(Pr2[thr][k]-Pr2[bot][k]-f3*(Pr2[top][k]-Pr2[bot][k]));
			}
		}
	}
}
double evpot(double R2,double z){//returns Phi from phil
	double pol[npoly],phip[npoly];
	double c=z/sqrt(R2+z*z), r=sqrt(R2+z*z);
	intpo2(r,phip,npoly); evenlegend(pol,c,npoly);
	//if(z==0. && r>1e-6 ) printf("### R=%f phip[0]=%f phip[1]=%f phip[2]=%f\n",r,phip[0],phip[1],phip[2]);
	double phi=phip[0];
	for(int np=1; np<npoly; np++) phi+=phip[np]*pol[np];
	return phi;
}
double dPhir(double r,double theta,double *d2Phidr2){
	double pol[npoly],phip[npoly],dphip[npoly],d2phip[npoly];
	double grad1=0,grad2=0;
	evenlegend(pol,cos(theta),npoly);
	intpo2(r,phip,dphip,d2phip,npoly);
	for(int np=0; np<npoly; np++){
		grad1+=dphip[np]*pol[np];
		grad2+=d2phip[np]*pol[np];
	}
	(*d2Phidr2)=grad2; return grad1;
}
double dPhitheta(double r, double c) {
    //returns the first derivative of the whole potential wrt spherical theta
	double dphitheta = 0;
	double dpol[npoly], phinterp[npoly];
	intpo2(r, phinterp, npoly);
	for(int k=0; k<npoly; k++) dpol[k] = dlegend(c, 2*k); //fills the array dpol with derivatives of the even legendre polys
	for(int j=0; j<npoly; j++) dphitheta += phinterp[j]*dpol[j]; //this is the derivative wrt cos(theta), need to then multiply by -sin(theta)
	dphitheta *= -sqrt(1 - c*c); //im taking the positive root here - any issues with this?
	return dphitheta;
}
double d2Phitheta(double r, double c) {
    //returns the second derivative of the whole potential wrt spherical theta
	double d2phitheta = 0;
	double d2pol[npoly], phinterp[npoly];
	intpo2(r, phinterp, npoly);
	for(int k=0; k<npoly; k++) d2pol[k] = d2legend(c, 2*k);
	for(int j=0; j<npoly; j++) d2phitheta += phinterp[j]*d2pol[j]; //2nd derivative wrt cos(theta)
	d2phitheta *= (1 - c*c);
	d2phitheta += ( (dPhitheta(r, c)*c) / sqrt(1 - c*c) ); //add the first derivative term from the product rule and ensure the correct factor
	return d2phitheta;
}
double phileg(double R,double z){
	if(R==0 && z==0)  R=ar[0];
	return evpot(R*R,z);
}
void get_Ylm(double **sig,double **sigl,int nr,int ngauss,int npoly){
	//calculates Legendre poly coeffs sigl of data in sig
	double ci[ngauss],wi[ngauss],pol[npoly],**poly;
	poly=dmatrix(npoly,ngauss);
	gauleg(0,1,ci,wi,ngauss);
	for(int i=0; i<ngauss; i++){
		evenlegend(pol,ci[i],npoly);
		for(int np=0; np<npoly; np++) poly[np][i]=2*pol[np]*wi[i];
	}
	for(int n=1; n<nr; n++){
		for(int np=0; np<npoly; np++){
			sigl[n][np]=0;
			for(int i=0; i<ngauss; i++)
				sigl[n][np]+=poly[np][i]*sig[n][i];
		}
	}
	sigl[0][0]=2*sig[0][0];
	for(int i=1; i<npoly; i++) sigl[0][i]=0;
	delmatrix(poly,npoly);
}
/*
 *  integrands I_int and I_ext for potential computation
 *  they use linear interpolation on rhl grid
 */
double I_int_f(double x, void *par){
	int l = *(int *) par;

	double rl[nr];
	for (int i=0; i<nr; i++) rl[i]=rhl[i][l/2];

	if (x<ar[nr-1]) return pow(x,l+2)*linterp(ar,rl,nr,x);
	else return 0.;
}
double I_ext_f(double x, void *par){
	int l = *(int *) par;

	double rl[nr];
	for (int i=0; i<nr; i++) rl[i]=rhl[i][l/2];

	if (x<ar[nr-1]) return pow(x,-l+1)*linterp(ar,rl,nr,x);
	else return MAX(0.,pow(x,-l+1)*( rl[nr-1]+(rl[nr-1]-rl[nr-2])/(ar[nr-1]-ar[nr-2])*(x-ar[nr-1]) ));
}
void potent(char *fname,double (*dens)(double,double),int rd,int prnt){
// Puts into rhl and phil  Legendre poly coefficients of density & potential associated
// with mass density dens; if rd==0 compute rhl, otherwise read it
	double **poly,**rho,**I_int,**I_ext;
	double pol[npoly];
	poly=dmatrix(npoly,ngauss); rho=dmatrix(nr,ngauss);
	I_int=dmatrix(nr,npoly); I_ext=dmatrix(nr,npoly);
	double ci[ngauss],si[ngauss],wi[ngauss];
	gauleg(0,1,ci,wi,ngauss);
// npoly coeffs for Gauss integration over angles ----------------------------
	for(int i=0; i<ngauss; i++){
		si[i]=sqrt(1-ci[i]*ci[i]);
		evenlegend(pol,ci[i],npoly);
		for(int np=0; np<npoly; np++) poly[np][i]=2*pol[np]*wi[i];
	}
	FILE *ofile;
	if(rd==0){//compute, don't read data
		ofile=fopen(fname,"w");
		if(ofile==NULL) printf("cannot open file %s\n",fname);

		/**************************************************************
		 *  the radial grid starts at TINY now..
		 */
		rho[0][0]=(*dens)(ar[0],0);
		printf("rho0 %g ",rho[0][0]);

		/* smooth truncation of density */
		for(int i=0; i<ngauss; i++) rho[nr-1][i]=(*dens)(ar[nr-1]*si[i],ar[nr-1]*ci[i]);

		for(int j=1; j<ngauss; j++) rho[0][j]=rho[0][0];
		for(int n=1; n<nr; n++){//Get density @ grid pts
			if(prnt==1) printf("%d ",n);
			double r=ar[n];
            for(int np=0; np<npoly; np++) rhl[n][np]=0;
#pragma omp parallel for
			for(int i=0; i<ngauss; i++){
				rho[n][i]=(*dens)(r*si[i],r*ci[i]) - rho[nr-1][i];	
            }
			for(int np=0; np<npoly; np++)
				for(int i=0; i<ngauss; i++){
					rhl[n][np]+=poly[np][i]*rho[n][i];
				}
		}
		rhl[0][0]=2*rho[0][0];
		for(int i=1;i<npoly;i++) rhl[0][i]=0;
		if(prnt==1) printf("\n");
		fprintf(ofile,"%d %d %d\n",nr,npoly,ngauss);
		compress(ofile,ar,nr); compres2(ofile,rhl,nr,npoly);
		compres2(ofile,rho,nr,ngauss);
	}else{
		ofile=fopen(fname,"r");
		if(ofile==NULL){
			printf("cannot open file %s\n",fname); exit(0);
		}
		if (fscanf(ofile,"%d %d %d",&nr,&npoly,&ngauss)==0) printf("[WARNING] read of %s not successful..\n",fname);
		get(ofile,ar,nr); get2(ofile,rhl,nr,npoly); get2(ofile,rho,nr,ngauss);
	}
	fclose(ofile);
// rhl[n][i] now contains l=2*i legendre coefficient of angular
// distribution of density on shell radius ar[n]
	int nc=nr/12;// change integration technique at this grid pt
	for(int i=0;i<npoly;i++){
		I_int[0][i]=0; I_ext[0][i]=0;//integrals from zero to zero
	}

#if defined(HERNQUIST) || defined(NFW)
	/*
	 * 	Corrections for singular models (e.g. Hernquist)
	 */
	for(int np=0; np<npoly; np++){
		int l=2*np;
		double r=ar[0],a,p1=0;
		/*
		 * get rid of quadratic terms in l=0
		 *
		 *	if(np==0){
		 *		a=rhl[0][np];
		 *		I_int[1][np]=-.5*(a*ar[0])*(pow(ar[1],2)-pow(ar[0],2));		// before: -.2*(a-rhl[1][0])*pow(ar[1],3);
		 *		I_ext[1][np]=-(a*ar[0])*(ar[1]-ar[0]);		// before: -.25*(a-rhl[1][0])*pow(ar[1],2);
		 *	} else {//assume const density inside ar[1] and rhl_l~r^l
		 *		a=rhl[1][np]/pow(ar[1],l-1);			// singular version, before:rhl[1][np]/pow(ar[1],l)
		 *	}
		*/

		/*
		 *  replacing these integrals with their singular r^-1 versions
		 */
		p1=pow(ar[0],2*l+2); a=rhl[0][np]/pow(ar[0],l-1);
		for(int n=1; n<nc; n++){//leading terms in all integrals
			r=ar[n]; double p2=pow(r,2*l+2);
			double b=rhl[n][np]/pow(r,l-1), arg=.5*(a+b);
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+2);
			I_ext[n][np]=I_ext[n-1][np]+arg*(r-ar[n-1]);
			a=b; p1=p2;
		}
		a*=pow(r,l); p1=pow(r,l+3);
		for(int n=nc; n<nr; n++){
			r=ar[n]; double p2=pow(r,l+3);
			double b=rhl[n][np], arg=.5*(a+b); 
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+3);
			if(l!=2) I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,2-l)-pow(ar[n-1],2-l))/(double)(2-l);
			else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar[n-1]));
			a=b; p1=p2;
		}
	}
#elif defined ISOTHERMAL
	/*
	 * 	Corrections for singular isothermal models

	for(int np=0; np<npoly; np++){
		int l=2*np;
		double r=ar[0],a,p1=0;

		//grid does not start at 0
		I_int[0][np]=.5*rhl[0][np]/pow(ar[0],l-2)*pow(ar[0],2*l+1)/(double)(2*l+1);

		//  replacing these integrals with their singular r^-2 versions
		p1=pow(ar[0],2*l+1); a=rhl[0][np]/pow(ar[0],l-2);
		for(int n=1; n<nr; n++){//leading terms in all integrals
			r=ar[n]; double p2=pow(r,2*l+1);
			double b=rhl[n][np]/pow(r,l-2), arg=.5*(a+b);
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+1);
			I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar[n-1]));
			a=b; p1=p2;
		}
	}
	*/

	/*
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	for(int np=0; np<npoly; np++){
			int l=2*np;

			gsl_function f; double err;

			for(int n=0; n<nr; n++){
				f.function = &I_int_f;
				f.params = &l;

				gsl_integration_qags(&f,0.,ar[n],1e-5,1e-4,1000,w,&I_int[n][np],&err);

				f.function = &I_ext_f;
				f.params = &l;
				gsl_integration_qagiu(&f,ar[n],1e-5,1e-4,1000,w,&I_ext[nr-1-n][np],&err);
			}
		}
	gsl_integration_workspace_free (w);
	*/


	for(int np=0; np<npoly; np++){
			int l=2*np;
			double r,a,p1=0;
			if(np==0){//corrections for quadratic term in a at a<a[1]
				a=rhl[0][np];
				I_int[1][np]=-.2*(a-rhl[1][0])*pow(ar[1],3);
				I_ext[1][np]=-.25*(a-rhl[1][0])*pow(ar[1],2);
			} else {//assume const density inside ar[1] and rhl_l~r^l
				a=rhl[1][np]/pow(ar[1],l);
			}
			for(int n=1; n<nc; n++){//leading terms in all integrals
				r=ar[n]; double p2=pow(r,2*l+3);
				double b=rhl[n][np]/pow(r,l), arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+3);
				I_ext[n][np]=I_ext[n-1][np]+.5*arg*(r*r-ar[n-1]*ar[n-1]);
				a=b; p1=p2;
			}
			// replaced for r^-2 at r --> + infty
			a=rhl[nc-1][np]*pow(ar[nc-1],2); p1=pow(ar[nc-1],l+1);
			for(int n=nc; n<nr; n++){
				r=ar[n]; double p2=pow(r,l+1);
				double b=rhl[n][np]*pow(r,2), arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+1);
				if(l!=0)I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,-l)-pow(ar[n-1],-l))/(double)(-l);
				else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar[n-1]));
				a=b; p1=p2;
			}
	}

#elif defined ISOCHRONE
	for(int np=0; np<npoly; np++){
			int l=2*np;
			double r,a,p1=0;
			if(np==0){//corrections for quadratic term in a at a<a[1]
				a=rhl[0][np];
				I_int[1][np]=-.2*(a-rhl[1][0])*pow(ar[1],3);
				I_ext[1][np]=-.25*(a-rhl[1][0])*pow(ar[1],2);
			} else {//assume const density inside ar[1] and rhl_l~r^l
				a=rhl[1][np]/pow(ar[1],l);
			}
			for(int n=1; n<nc; n++){//leading terms in all integrals
				r=ar[n]; double p2=pow(r,2*l+3);
				double b=rhl[n][np]/pow(r,l), arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+3);
				I_ext[n][np]=I_ext[n-1][np]+.5*arg*(r*r-ar[n-1]*ar[n-1]);
				a=b; p1=p2;
			}
			a*=pow(r,l); p1=pow(r,l+3);
			for(int n=nc; n<nr; n++){
				r=ar[n]; double p2=pow(r,l+3);
				double b=rhl[n][np], arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+3);
				if(l!=2) I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,2-l)-pow(ar[n-1],2-l))/(double)(2-l);
				else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar[n-1]));
				a=b; p1=p2;
			}
		}
#endif
	for(int np=0; np<npoly; np++){
		int l=2*np;
		if(l==0){
			//phil[0][np]=-TPI*G*(I_ext[nr-1][np]-I_ext[0][np]);
			phil[0][np]=-TPI*G*(I_int[0][np]/pow(ar[0],l+1)+
						(I_ext[nr-1][np]-I_ext[0][np])*pow(ar[0],l));
			Pr2[0][np]=TPI*G*rhl[0][np];
		}else{
			phil[0][np]=0;
			if(l==2) Pr2[0][np]=-TPI*G*((l-1)*l*I_ext[nr-1][np]-(2*l+1)*rhl[0][np]);
			else Pr2[0][np]=TPI*G*(2*l+1)*rhl[0][np];
		}
		Pr[0][np]=0;
		for(int n=1; n<nr; n++){
			phil[n][np]=-TPI*G*(I_int[n][np]/pow(ar[n],l+1)+
					    (I_ext[nr-1][np]-I_ext[n][np])*pow(ar[n],l));
			Pr[n][np]=-TPI*G*(-(l+1)*I_int[n][np]/pow(ar[n],l+2)+
					  l*(I_ext[nr-1][np]-I_ext[n][np])*pow(ar[n],l-1));
			Pr2[n][np]=-TPI*G*((l+2)*(l+1)*I_int[n][np]/pow(ar[n],l+3)+
					   (l-1)*l*(I_ext[nr-1][np]-I_ext[n][np])*pow(ar[n],l-2)
					   -(2*l+1)*rhl[n][np]);
		}
	}
	/* output of rho-Phi profiles */
	char times[8],orho[30]; sprintf(times,"_%d.out",iter);
	strcpy(orho,"rhof"); strcat(orho,times);
	FILE* rhof=fopen(orho,"w");
	for (int i=0;i<nr;i++) {
		fprintf(rhof,"%e %e %e %e %e %e\n",ar[i],rho[i][0],evpot(ar[i]*ar[i],0),phi_isoth(ar[i]),rho_isoth(ar[i]),phil[i][0]);
	}

	// debug purposes
	for (int n=0;n<nr;n++) printf(">> phil[%d][0]=%f rhl[%d][0]=%f rho[%d][0]=%f\n",n,phil[n][0],n,rhl[n][0],n,rho[n][0]);
	fclose(rhof); iter++;	delmatrix(poly,npoly); delmatrix(rho,nr);
	delmatrix(I_int,nr); delmatrix(I_ext,nr);
}
void potent5(char *fname,double (*dens)(double,double,double*,double*,double*,double*,double*,double),int rd,int prnt){
// Puts into rhl and phil Legendre poly coefficients of density & potential associated
// with mass density dens; if rd==0 compute rhl, otherwise read it
	double **poly,**rho,**Vbar,**sigR,**sigp,**sigz,**sigRz,**I_int,**I_ext;
	double pol[npoly];
	poly=dmatrix(npoly,ngauss);
	rhl=dmatrix(nr,npoly); rho=dmatrix(nr,ngauss); Vbar=dmatrix(nr,ngauss);
	sigR=dmatrix(nr,ngauss); sigp=dmatrix(nr,ngauss);
	sigz=dmatrix(nr,ngauss); sigRz=dmatrix(nr,ngauss);
	I_int=dmatrix(nr,npoly); I_ext=dmatrix(nr,npoly);
	double ci[ngauss],si[ngauss],wi[ngauss];
	gauleg(0,1,ci,wi,ngauss);
	//for(int i=0;i<ngauss;i++) printf("%f ",2*wi[i]);
// npoly coeffs for Gauss integration over angles ----------------------------
	for(int i=0; i<ngauss; i++){
		si[i]=sqrt(1-ci[i]*ci[i]);
		evenlegend(pol,ci[i],npoly);
		for(int np=0; np<npoly; np++) poly[np][i]=2*pol[np]*wi[i];
	}
	FILE *ofile;
	if(rd==0){//compute, don't read data
		ofile=fopen(fname,"w");
		printf("--starting rho evaluation..\n");
		if(ofile==NULL) printf("cannot open file %s\n",fname);

		/*	******************************************************
		 *  Performing Backwards integration in radial grid
		 *
		 *  Update 19/02: actually there's no difference if using
		 *  sinh-integration or standard and neither if setting
		 *  Vlim to a value <V_escape.
		 *  This is now set as before, only the backward integration
		 *  is kept.
		 */

		/* smooth truncation of density */
		for(int i=0; i<ngauss; i++)
			rho[nr-1][i]=(*dens)(ar[nr-1]*si[i],ar[nr-1]*ci[i],Vbar[nr-1]+i,
				sigR[nr-1]+i,sigp[nr-1]+i,sigz[nr-1]+i,sigRz[nr-1]+i,
				sqrt(-2*(Phi(ar[nr-1]*si[i],ar[nr-1]*ci[i])-Phi(100,100))));

		/* initialize bar */
		ProgressBar bar(60);
		bar.init(nr);
		int nn=0;

		for(int n=nr-1; n>0; n--){//Get density @ grid pts
			/* update bar*/
			bar.update(nn+1); nn++;

			//if(prnt==1) printf("%d ",n);
			double Vscale;
			double r=ar[n];
			for(int np=0; np<npoly; np++) rhl[n][np]=0;
#pragma omp parallel for
			for(int i=0; i<ngauss; i++){
				// OLDER THAN 19/02: control on the velocity integration-limits

				   // for the innermost grid points I set the integration limit to be 2*sigma[k-1]
				 	if (n<nr/10){
				 		Vscale=2.*sqrt(pow(*(sigR[n+1]+i),2)+pow(*(sigp[n+1]+i),2)+pow(*(sigz[n+1]+i),2));
				 		//Vscale=sqrt(-2*(Phi(r*si[i],r*ci[i])-Phi(100,100)));
				 		rho[n][i]=(*dens)(r*si[i],r*ci[i],Vbar[n]+i,sigR[n]+i,sigp[n]+i,sigz[n]+i,sigRz[n]+i,Vscale) - rho[nr-1][i];
				 	}
				 	else {
				 		Vscale=sqrt(-2*(Phi(r*si[i],r*ci[i])-Phi(100,100)));
				 		rho[n][i]=(*dens)(r*si[i],r*ci[i],Vbar[n]+i,sigR[n]+i,sigp[n]+i,sigz[n]+i,sigRz[n]+i,Vscale) - rho[nr-1][i];
				 	}


				//rho[n][i]=(*dens)(r*si[i],r*ci[i],Vbar[n]+i,sigR[n]+i,sigp[n]+i,sigz[n]+i,sigRz[n]+i);
				if(isnan(rho[n][i])){
					rho[n][i]=(*dens)(r*si[i],r*ci[i],Vbar[n]+i,sigR[n]+i,sigp[n]+i,sigz[n]+i,sigRz[n]+i,Vscale);
					printf("in potent5: %d %d %g\n",n,i,rho[n][i]);
					exit(0);
				}
			}
			//printf("%f %f %f %f\n",ar[n],Vbar[n][0],Vbar[n][1],Vbar[n][2]);
		}
		/* finalize bar */
		bar.fillSpace("..done potential computation!!\n\n");

		/**************************************************************
		 *  the radial grid starts at TINY now..
		 */
		double Vscale=2.*sqrt(pow(*(sigR[1]),2)+pow(*(sigp[1]),2)+pow(*(sigz[1]),2));
		rho[0][0]=(*dens)(ar[0],0,Vbar[0],sigR[0],sigp[0],sigz[0],sigRz[0],Vscale);
		//sigp[0][0]=sqrt(sigp[0][0]+pow(Vbar,2));

		/*
		 * 03/09/14: linear interpolation in the first grid point!
		 * 			 it was causing troubles in the potential computation
		 * 			 since the density was extremely high, implying very high Phi
		 */
		rho[0][0] = (rho[2][0]-rho[1][0])/(ar[2]-ar[1])*(ar[0]-ar[1])+rho[1][0];
		printf("rho0 %g %f %f\n",rho[0][0],sigR[0][0],Vbar[0][0]);
		for(int j=1; j<ngauss; j++){
			rho[0][j]=rho[0][0]; sigR[0][j]=sigR[0][0];
			sigp[0][j]=sigp[0][0]; sigz[0][j]=sigz[0][0];
			sigRz[0][j]=sigRz[0][0]; Vbar[0][j]=0;
		}

		/* Computing rhl */
		for (int n=0; n<nr; n++)
			for(int np=0; np<npoly; np++)
				for(int i=0; i<ngauss; i++)
					rhl[n][np]+=poly[np][i]*rho[n][i];


		rhl[0][0]=2*rho[0][0];
		/********************************************************************************
		 * 	BUG in the code! was "i<ngauss", now set to "i<npoly"
		 */
		for(int i=1;i<npoly;i++) rhl[0][i]=0;
		if(prnt==1) printf("\n");
		fprintf(ofile,"%d %d %d\n",nr,npoly,ngauss);
		compress(ofile,ar,nr); compres2(ofile,rhl,nr,npoly);
		compres2(ofile,rho,nr,ngauss); compres2(ofile,sigR,nr,ngauss);
		compres2(ofile,sigp,nr,ngauss); compres2(ofile,sigz,nr,ngauss);
		compres2(ofile,sigRz,nr,ngauss); compres2(ofile,Vbar,nr,ngauss);
	}else{
		ofile=fopen(fname,"r");
		if(ofile==NULL){
			printf("cannot open file %s\n",fname); exit(0);
		}
		if (fscanf(ofile,"%d %d %d",&nr,&npoly,&ngauss)==0) printf("[WARNING] read of %s not successful..\n",fname);
		get(ofile,ar,nr); get2(ofile,rhl,nr,npoly);
	}
	fclose(ofile);
// rhl[n][i] now contains l=2*i legendre coefficient of angular
// distribution of density on shell radius ar[n]
	int nc=nr/12;// change integration technique at this grid pt
	for(int i=0;i<npoly;i++){
		I_int[0][i]=0; I_ext[0][i]=0;//integrals from zero to zero
	}

#if defined(HERNQUIST) || defined(NFW)
	/*
	 * 	Corrections for singular models (e.g. Hernquist)
	 */
	for(int np=0; np<npoly; np++){
		int l=2*np;
		double r=ar[0],a,p1=0;
		/*
		 * get rid of quadratic terms in l=0
		 *
		 *	if(np==0){
		 *		a=rhl[0][np];
		 *		I_int[1][np]=-.5*(a*ar[0])*(pow(ar[1],2)-pow(ar[0],2));		// before: -.2*(a-rhl[1][0])*pow(ar[1],3);
		 *		I_ext[1][np]=-(a*ar[0])*(ar[1]-ar[0]);		// before: -.25*(a-rhl[1][0])*pow(ar[1],2);
		 *	} else {//assume const density inside ar[1] and rhl_l~r^l
		 *		a=rhl[1][np]/pow(ar[1],l-1);			// singular version, before:rhl[1][np]/pow(ar[1],l)
		 *	}
		*/

		/*
		 *  replacing these integrals with their singular r^-1 versions
		 */
		p1=pow(ar[0],2*l+2); a=rhl[0][np]/pow(ar[0],l-1);
		for(int n=1; n<nc; n++){//leading terms in all integrals
			r=ar[n]; double p2=pow(r,2*l+2);
			double b=rhl[n][np]/pow(r,l-1), arg=.5*(a+b);
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+2);
			I_ext[n][np]=I_ext[n-1][np]+arg*(r-ar[n-1]);
			a=b; p1=p2;
		}
		a*=pow(r,l); p1=pow(r,l+3);
		for(int n=nc; n<nr; n++){
			r=ar[n]; double p2=pow(r,l+3);
			double b=rhl[n][np], arg=.5*(a+b);
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+3);
			if(l!=2) I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,2-l)-pow(ar[n-1],2-l))/(double)(2-l);
			else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar[n-1]));
			a=b; p1=p2;
		}
	}
#elif defined ISOTHERMAL
	/*
	 * 	Corrections for singular isothermal models

	for(int np=0; np<npoly; np++){
		int l=2*np;
		double r=ar[0],a,p1=0;

		//grid does not start at 0
		I_int[0][np]=.5*rhl[0][np]/pow(ar[0],l-2)*pow(ar[0],2*l+1)/(double)(2*l+1);

		//  replacing these integrals with their singular r^-2 versions
		p1=pow(ar[0],2*l+1); a=rhl[0][np]/pow(ar[0],l-2);
		for(int n=1; n<nr; n++){//leading terms in all integrals
			r=ar[n]; double p2=pow(r,2*l+1);
			double b=rhl[n][np]/pow(r,l-2), arg=.5*(a+b);
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+1);
			I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar[n-1]));
			a=b; p1=p2;
		}
	}
	*/
	for(int np=0; np<npoly; np++){
			int l=2*np;
			double r,a,p1=0;
			if(np==0){//corrections for quadratic term in a at a<a[1]
				a=rhl[0][np];
				I_int[1][np]=-.2*(a-rhl[1][0])*pow(ar[1],3);
				I_ext[1][np]=-.25*(a-rhl[1][0])*pow(ar[1],2);
			} else {//assume const density inside ar[1] and rhl_l~r^l
				a=rhl[1][np]/pow(ar[1],l);
			}
			for(int n=1; n<nc; n++){//leading terms in all integrals
				r=ar[n]; double p2=pow(r,2*l+3);
				double b=rhl[n][np]/pow(r,l), arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+3);
				I_ext[n][np]=I_ext[n-1][np]+.5*arg*(r*r-ar[n-1]*ar[n-1]);
				a=b; p1=p2;
			}/*
			a*=pow(r,l); p1=pow(r,l+3);
			for(int n=nc; n<nr; n++){
				r=ar[n]; double p2=pow(r,l+3);
				double b=rhl[n][np], arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+3);
				if(l!=2) I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,2-l)-pow(ar[n-1],2-l))/(double)(2-l);
				else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar[n-1]));
				a=b; p1=p2;
			}
			*/
			// replaced for r^-2 at r --> + infty
			a=rhl[nc-1][np]*pow(ar[nc-1],2); p1=pow(ar[nc-1],l+1);
			for(int n=nc; n<nr; n++){
				r=ar[n]; double p2=pow(r,l+1);
				double b=rhl[n][np]*pow(r,2), arg=.5*(a+b);
				I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+1);
				if(l!=0)I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,-l)-pow(ar[n-1],-l))/(double)(-l);
				else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar[n-1]));
				a=b; p1=p2;
			}
	}
#elif defined ISOCHRONE
	for(int np=0; np<npoly; np++){
		int l=2*np;
		double r,a,p1=0;
		if(np==0){//corrections for quadratic term in a at a<a[1]
			a=rhl[0][np];
			I_int[1][np]=-.2*(a-rhl[1][0])*pow(ar[1],3);
			I_ext[1][np]=-.25*(a-rhl[1][0])*pow(ar[1],2);
		} else {//assume const density inside ar[1] and rhl_l~r^l
			a=rhl[1][np]/pow(ar[1],l);
		}
		for(int n=1; n<nc; n++){//leading terms in all integrals
			r=ar[n]; double p2=pow(r,2*l+3);
			double b=rhl[n][np]/pow(r,l), arg=.5*(a+b); 
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(2*l+3);
			I_ext[n][np]=I_ext[n-1][np]+.5*arg*(r*r-ar[n-1]*ar[n-1]);
			a=b; p1=p2;
		}
		a*=pow(r,l); p1=pow(r,l+3);
		for(int n=nc; n<nr; n++){
			r=ar[n]; double p2=pow(r,l+3);
			double b=rhl[n][np], arg=.5*(a+b); 
			I_int[n][np]=I_int[n-1][np]+arg*(p2-p1)/(double)(l+3);
			if(l!=2) I_ext[n][np]=I_ext[n-1][np]+arg*(pow(r,2-l)-pow(ar[n-1],2-l))/(double)(2-l);
			else I_ext[n][np]=I_ext[n-1][np]+arg*(log(r)-log(ar[n-1]));
			a=b; p1=p2;
		}
	}
#endif
	for(int np=0; np<npoly; np++){
		int l=2*np;
		if(l==0){
			phil[0][np]=-TPI*G*(I_ext[nr-1][np]-I_ext[0][np]);
			Pr2[0][np]=TPI*G*rhl[0][np];
		}else{
			phil[0][np]=0;
			if(l==2) Pr2[0][np]=-TPI*G*((l-1)*l*I_ext[nr-1][np]-(2*l+1)*rhl[0][np]);
			else Pr2[0][np]=TPI*G*(2*l+1)*rhl[0][np];
		}
		Pr[0][np]=0;
		for(int n=1; n<nr; n++){
			phil[n][np]=-TPI*G*(I_int[n][np]/pow(ar[n],l+1)+
					    (I_ext[nr-1][np]-I_ext[n][np])*pow(ar[n],l));
			Pr[n][np]=-TPI*G*(-(l+1)*I_int[n][np]/pow(ar[n],l+2)+
					  l*(I_ext[nr-1][np]-I_ext[n][np])*pow(ar[n],l-1));
			Pr2[n][np]=-TPI*G*((l+2)*(l+1)*I_int[n][np]/pow(ar[n],l+3)+
					   (l-1)*l*(I_ext[nr-1][np]-I_ext[n][np])*pow(ar[n],l-2)
					   -(2*l+1)*rhl[n][np]);
		}
	}
	/* output of rho-Phi profiles */
	char times[8],orho[30]; sprintf(times,"_%d.out",iter);
	strcpy(orho,"rhof"); strcat(orho,times);
	FILE* rhof=fopen(orho,"w");
	for (int i=0;i<nr;i++) {
		fprintf(rhof,"%e %e %e %e %e %e\n",ar[i],rho[i][0],evpot(ar[i]*ar[i],0),phi_isoth(ar[i]),rho_isoth(ar[i]),phil[i][0]);
	}
	fclose(rhof); iter++;

	// debug purposes
	for (int n=0;n<nr;n++) printf(">> phil[%d][0]=%f rhl[%d][0]=%f rho[%d][0]=%f\n",n,phil[n][0],n,rhl[n][0],n,rho[n][0]);
	// check forces: -G*M/r^2=Pr(r) for very large r (Keplerian regime)
	//for (int n=0;n<nr;n++) printf(">> -1/ar^2[%d]=%f Pr[%d][0]=%f Pr2[%d][0]=%f\n",n,-1./pow(ar[n],2),n,Pr[n][0],n,Pr2[n][0]);
	delmatrix(poly,npoly); delmatrix(rho,nr);
	delmatrix(sigR,nr); delmatrix(sigp,nr); delmatrix(sigz,nr);
	delmatrix(sigRz,nr);
	delmatrix(I_int,nr); delmatrix(I_ext,nr);
}
