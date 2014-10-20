#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uv_orb3.h"
#include "tables3.h"
#include "press.h"
#include "ELz_grid.hh"
#include "eddington.hh"
#include "moments.h"
#include "leg_pot2.h"
#include "ini_potLEG.h"
#include "oct_int.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_interp.h"

#define TPI 6.283185307179586
#define SQRT2 1.4142135623730951
#define SMALLEST 10		// see tables3.cpp
#define DFGRID 200

static const double GM=1,a=1;
static double aphi,az;
/* array of DF discretized on Egrid */
static double *df_coarse,*E_grid;
double *EofJ=NULL,*LzofJ;
static int gR,count,is_EofJ=0,is_DFgrid=0,ndf=0;
static FILE* fdf,*fdf2;

double Jr_iso(double E, double L){
	return GM/sqrt(-2*E)-.5*(L+sqrt(L*L+4*GM*a));
	double Et=-2*E*a/GM,l=L/(2*sqrt(GM*a));
	return -sqrt(GM*a)*(pow(Et,-1/2)-(sqrt(l*l+1)+l));
}
double get_j(double E){//barycentre of const-H surface of isochrone
	return (4*GM/sqrt(-2*E)-sqrt(4*pow(GM,2)/(-2*E)+12*GM*a))/6;
}
double H_iso(double Jr,double Lz,double Jz){//Hamiltonian of ischrone
	double L=fabs(Lz)+Jz;
	return -.5*pow(GM,2)/pow(Jr+.5*(L+sqrt(L*L+4*GM*a)),2);
}
void Omega_iso(double Jr,double L,double *Omegar,double *OmegaL){
	*Omegar=pow(GM,2)/pow(Jr+.5*(L+sqrt(L*L+4*GM*a)),3);
	*OmegaL=.5*(1+L/sqrt(L*L+4*GM))*(*Omegar);
}	
void setdf(double aphi_in,double az_in){
	aphi=aphi_in; az=az_in;
}
double df_iso(double Et){//isochr DF(dimless E) argument curly Etilde shld be positive
	if(Et<0 || Et>.5) return 0;
	double fac=pow(TPI,3)*pow(2*GM*a,1.5)/2;
	double Et5=sqrt(Et);
	return Et5/fac/pow(2.*(1.-Et),4)*(27+Et*(-66+Et*(320+Et*(-240+64*Et)))
				    +3*asin(Et5)/sqrt(Et*(1-Et))*(-9+Et*(28+16*Et)));
}
double df_Hern(double E){
	if(E<0 || E>1.) return 0.;
	double fac=SQRT2*pow(TPI,3)*pow(GM*a,1.5);
	return sqrt(E)/pow(1-E,2)/fac*((1-2.*E)*(8.*E*E-8.*E-3.)+3.*asin(sqrt(E))/sqrt(E*(1-E)));
}
/*
 *  Discretized DF of the isochrone model. E_grid from class ELz_grid
 */
void build_df_grid(void){

	if (is_DFgrid==0){
		ELz_grid grid;
		grid.make_NRgrids();
		gR=grid.gR;

		double *rg=dmatrix(DFGRID);
		df_coarse=dmatrix(DFGRID); E_grid=dmatrix(DFGRID);

		for (int i=0; i<DFGRID; i++)
			rg[i]=pow(10.,(log10(ar[0]) + (log10(ar[nr-1])-log10(ar[0]))*i/(double)(DFGRID-1) ));

		Edd_DF edd_df(rg,DFGRID);
		for (int i=0; i<DFGRID; i++) {
				df_coarse[i]= edd_df.DF[i]; 		// I calculate the discretized DF on
				E_grid[i]=edd_df.psi_r[i];			// -E*a/GM, while the grid is in E
		}

		is_DFgrid=1;
	}
}

void read_EofJ(void){
	int c;
	rewind(fdf);

	/* count lines */
	while ( (c=fgetc(fdf)) != EOF ) {
        if ( c == '\n' )
            ndf++;
    }
	rewind(fdf);

	double dummy;
	EofJ=dmatrix(ndf); LzofJ=dmatrix(ndf);
	for (int i=0;i<ndf;i++) {
		if (fscanf(fdf,"%lf %lf %lf %lf %lf\n",&LzofJ[i],
	   			  &dummy,&EofJ[i],&dummy,&dummy)==0)
		    	   printf("\n[WARNING]: cannot read df.dat\n");
		EofJ[i]=-EofJ[i];
	}

	sort(ndf,LzofJ); sort(ndf,EofJ);
	is_EofJ=1;
}

void build_EofJ(void){

	if (is_EofJ==0){
		fdf=fopen("df.dat","w+"); count=0;
		fdf2=fopen("en.dat","w+");

		/*
		 *  radial grid starts at tiny now..
		 */
		printf("..Tabulating Lz(J) and E(J) for spherical model, rho[0][0]=%f rho[%d][0]=%f\n",rho(ar[0],0),nr-1,rho(ar[nr-1],0));
		read_EofJ();

		//fclose(fdf);
		//remove("df.dat");
	}
}

/*
 * 	Distribution function for Isochrone potential as a function of the actions
 */
double df_I(double Jr,double Lz,double Jz,double R, double Etrue){
	double L=fabs(Lz)+Jz;

	/*
	 *	Isochrone case: use analytic formulae
	 *
	 *double H=-.5*pow(GM,2)/pow(Jr+.5*(L+sqrt(L*L+4*GM*a)),2);
	 *return df(-H*a/GM);
	 */

	double quadE(double Jr, double Lz, double Jz, double **surfLib, double * OmRat, int nsurfs, double Etrue);
	double H;
	if (is_EofJ==0){
		count++;

		ELz_grid grid;
		double Omegar,Omegaphi;
		Omega_iso(Jr,L,&Omegar,&Omegaphi);
		Omegaphi*=Lz/fabs(Lz);	// Omegaphi=sgn(Lz)*OmegaL
		double j=Omegar/Omegaphi*Jr+Lz+Lz/fabs(Lz)*Jz;
		//double j=Omegar/Omegaphi*Jr+fabs(Lz)+Jz;
		double Rc=grid.GetRc(j*j,R);
		H=.5*Rc*dPhiR(Rc)+PhiR(Rc);
		//if(count%1==0) fprintf(fdf,"%e %e %e %e %e %e %e\n",Lz,j,-H,.5*pow(GM,2)/pow(Jr+.5*(L+sqrt(L*L+4*GM*a)),2),
		//						linterp(E_grid,df_coarse,DFGRID,.5*pow(GM,2)/pow(Jr+.5*(L+sqrt(L*L+4*GM*a)),2)),
		//						df_iso(.5*pow(GM,2)/pow(Jr+.5*(L+sqrt(L*L+4*GM*a)),2)),linterp(E_grid,df_coarse,DFGRID,-H));

		if(count%1==0) fprintf(fdf,"%lf %lf %lf %lf %lf\n",Lz,j,-H,df_iso(-H),linterp(E_grid,df_coarse,DFGRID,-H));
	}
	else{
		//H=linterp(LzofJ,EofJ,ndf,Lz);
		double OmegaRatios[3]={0.,0.,0.},H = quadE(Jr,Lz,Jz,surf_LIB,OmegaRatios,nsurf,Etrue);
		return df_iso(-H*a/GM);
	}

	if (-H<0 || -H>.5) return 0.;

	return linterp(E_grid,df_coarse,DFGRID,-H);
}

/*
 * 	Distribution function for Hernquist potential as a function of the actions
 */
double df_H(double Jr,double Lz,double Jz,double *x, double *v){

	double H;
	if (is_EofJ==0){
		count++;

		double R=x[0],Lzh=R*v[1],Phigl=x[2];
		/* Hamiltonian: used first to compute Delta */
		H=.5*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2))+Phigl;

		double p[2]={v[0],v[2]};
		uv_orb uvorb(Deltafn(H),Lzh,Phigl,x,p);

		ELz_grid grid;
		uvorb.GetFreqs();
		double Omegar=uvorb.Omegau,Omegaphi=Lz/fabs(Lz)*uvorb.Omegav;
		double j=Omegar/Omegaphi*Jr+Lz+Lz/fabs(Lz)*Jz;
		//double j=Omegar/Omegaphi*Jr+fabs(Lz)+Jz;
		double Rc=grid.GetRc(j*j,R);
		H=.5*Rc*dPhiR(Rc)+PhiR(Rc);
		//if(count%1==0) fprintf(fdf,"%e %e %e %e %e %e %e\n",Lz,j,-H,.5*pow(GM,2)/pow(Jr+.5*(L+sqrt(L*L+4*GM*a)),2),
		//						linterp(E_grid,df_coarse,DFGRID,.5*pow(GM,2)/pow(Jr+.5*(L+sqrt(L*L+4*GM*a)),2)),
		//						df_iso(.5*pow(GM,2)/pow(Jr+.5*(L+sqrt(L*L+4*GM*a)),2)),linterp(E_grid,df_coarse,DFGRID,-H));

		if(count%1==0) fprintf(fdf,"%lf %lf %lf %lf %lf\n",Lz,j,-H,df_Hern(-H),linterp(E_grid,df_coarse,DFGRID,-H));
	}
	else{
		H=linterp(LzofJ,EofJ,ndf,Lz);
		//H-=phi_Hern(ar[nr-1])-evpot(ar[nr-1]*ar[nr-1],0.);
	}

	// Hernquist model: max E=1.
	if (-H<0 || -H>1.) return 0.;


	if (-H>E_grid[DFGRID-1]) {
		printf("-H=%g Egrid=%g\n",-H,E_grid[DFGRID-1]);
		double y2[DFGRID];
		spline(E_grid,df_coarse,DFGRID, 1., 0., y2);
		return splint(E_grid,df_coarse,y2,DFGRID,-H);
	}
	return linterp(E_grid,df_coarse,DFGRID,-H);
}

/*
 *  expdf_I: expdf for Isochrone potential
 */
double expdf_I(double j,double Jr,double Lz,double Jz,double R,
             double deltar,double deltaphi,double deltaz,double Etrue){
    double arg=+deltar*Jr+deltaphi*Lz+deltaz*Jz;
    //double delta=sqrt(deltar*deltar+deltaphi*deltaphi+deltaz*deltaz);
    return df_I(Jr,Lz,Jz,R,Etrue)*exp(-arg/j);
}
/*
 *  expdf_H: expdf for Hernquist potential
 */
double expdf_H(double j,double Jr,double Lz,double Jz,
			 double *x, double *v,
             double deltar,double deltaphi,double deltaz){
    double arg=+deltar*Jr+deltaphi*Lz+deltaz*Jz;
    //double delta=sqrt(deltar*deltar+deltaphi*deltaphi+deltaz*deltaz);
    return df_H(j,j,j,x,v)*exp(-arg/j);
}

double getJphi(double Jr, double Jz, double *a){
	return (sqrt((pow(a[8],2)-4*a[5]*a[6])*pow(Jz,2)+
			((2*a[8]*a[9]-4*a[6]*a[7])*Jr+2*a[3]*a[8]-4*a[2]*a[6])*Jz+
			(pow(a[9],2)-4*a[4]*a[6])*pow(Jr,2)+(2*a[3]*a[9]-4*a[1]*a[6])*Jr+4*a[0]*a[6]+pow(a[3],2))-
			a[8]*Jz-a[9]*Jr-a[3])/(2*a[6]);
}
double getJr(double Lz, double Jz, double *a){
	return (sqrt((pow(a[9],2)-4*a[4]*a[6])*pow(Lz,2)+
			((2*a[7]*a[9]-4*a[4]*a[8])*Jz+2*a[1]*a[9]-4*a[3]*a[4])*Lz+
			(pow(a[7],2)-4*a[4]*a[5])*pow(Jz,2)+(2*a[1]*a[7]-4*a[2]*a[4])*Jz+4*a[0]*a[4]+pow(a[1],2))-
			a[9]*Lz-a[7]*Jz-a[1])/(2*a[4]);
}
double getJz(double Jr, double Lz, double *a){
	return (sqrt((pow(a[8],2)-4*a[5]*a[6])*pow(Lz,2)+
			((2*a[7]*a[8]-4*a[5]*a[9])*Jr+2*a[2]*a[8]-4*a[3]*a[5])*Lz+
			(pow(a[7],2)-4*a[4]*a[5])*pow(Jr,2)+(2*a[2]*a[7]-4*a[1]*a[5])*Jr+4*a[0]*a[4]+pow(a[2],2))-
			a[8]*Lz-a[7]*Jr-a[2])/(2*a[5]);
}
double getJrPlane(double Lz, double Jz, double *a){
	return -(a[3]*Lz+a[2]*Jz-a[0])/a[1];
}

double normJacob(double Jr, double Jz, double Lz, double * a){
	return sqrt( pow(a[9]*Lz+a[7]*Jz+2.*a[4]*Jr+a[1],2)+pow(a[8]*Lz+2.*a[5]*Jz+a[7]*Jr+a[2],2)+pow(2.*a[6]*Lz+a[8]*Jz+a[9]*Jr+a[3],2) );
}
double ESurf(double Jr, double Jz, double Lz, double *a){
	return a[1]*Jr+a[2]*Jz+a[3]*Lz+a[4]*Jr*Jr+a[5]*Jz*Jz+a[6]*Lz*Lz+a[7]*Jr*Jz+a[8]*Jz*Lz+a[9]*Jr*Lz;
}
/*
 *  weighted average of the Energy corresponding to a given surface.
 *  The weights are the difference between the Etilde(J) and E(J).
 */
double ESurf(int surf, int nsurf, double Jr, double Jz, double Lz, double **surfLib){
	double dl= ( surf==0       ? 1e10 : fabs(ESurf(Jr,Jz,Lz,&surfLib[surf-1][0])-surfLib[surf-1][0]) ),
		   du= ( surf==nsurf-1 ? 1e10 : fabs(ESurf(Jr,Jz,Lz,&surfLib[surf+1][0])-surfLib[surf+1][0]) );

	double s2= MIN(dl,du); int surf2=(s2==dl ? surf-1 : surf+1);
	//double w1=1-fabs(ESurf(Jr,Jz,Lz,&surfLib[surf][0])-surfLib[surf][0]),w2=1-fabs(ESurf(Jr,Jz,Lz,&surfLib[surf2][0])-surfLib[surf2][0]);
	double w1=fabs(ESurf(Jr,Jz,Lz,&surfLib[surf][0])-surfLib[surf][0])   / normJacob(Jr,Jz,Lz,&surfLib[surf][0]),
		   w2=fabs(ESurf(Jr,Jz,Lz,&surfLib[surf2][0])-surfLib[surf2][0]) / normJacob(Jr,Jz,Lz,&surfLib[surf2][0]);

	return (w1*surfLib[surf][0]+
			w2*surfLib[surf2][0])/(w1+w2);
}

struct gslInterpFunctionParams { double *x,*y; int nn;};
double gslInterpFunction(double E, void * params){

	struct gslInterpFunctionParams * p = (gslInterpFunctionParams*) params;
	double *x=p->x, *y=p->y;
	int nn=p->nn;

	const gsl_interp_type *interpType;
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	interpType = gsl_interp_cspline;

	gsl_interp * interp = gsl_interp_alloc(interpType, nn);
	gsl_interp_init(interp, x, y, nn);

	double res = gsl_interp_eval_deriv (interp, x, y, E, acc);

	gsl_interp_free (interp);
	gsl_interp_accel_free (acc);

	return res;
}
double quadE(double Jr, double Lz, double Jz, double **surfLib, double * OmRat, int nsurfs, double Etrue){

	/*
	 * compute E from a_i coefficients of all quadratic surfaces
	 * and find closest to the real E of that surface.
	 */
	int surf=0,nn=0; double dd,diff=1e10;
	double EInterp[nsurfs],diffInterp[nsurfs];
	double E;
	for (int i=0; i<nsurfs; i++){

		dd=fabs(ESurf(Jr,Jz,Lz,&surfLib[i][0])-surfLib[i][0])/normJacob(Jr,Jz,Lz,&surfLib[i][0]);
		EInterp[nn]=-surfLib[i][0]; diffInterp[nn]=dd;
		if (i==0) nn++;
		else if (i>0 && EInterp[nn]>=EInterp[nn-1]) nn++;
		if ( dd<diff ) {
			surf=i; diff=dd;
		}
	}
	//if (diff/ESurf(Jr,Jz,Lz,&surfLib[surf][0])<-.1) return Etrue;

	//fprintf(fdf,"%e %e %e %e %e %e\n",diff,Etrue,ESurf(Jr,Jz,Lz,&surfLib[surf][0]),ESurf(surf,nsurfs,Jr,Jz,Lz,&surfLib[0]),Jr,getJrPlane(Lz,Jz,&surfLib[surf][0]));
	//fprintf(fdf,"%f %f %f\n",Etrue,ESurf(Jr,Jz,Lz,&surfLib[surf][0])-Etrue,ESurf(Jr,Jz,Lz,&surfLib[surf][0]));
	//printf("%d %f %f %f\n",surf,ESurf(surf,nsurfs,Jr,Jz,Lz,&surfLib[0]),Etrue,fabs(ESurf(surf,nsurfs,Jr,Jz,Lz,&surfLib[0])-Etrue));
	E=ESurf(surf,nsurfs,Jr,Jz,Lz,&surfLib[0]);

	OmRat[0]=surfLib[surf][2]/surfLib[surf][1];	// Omegaz   / Omegar
	OmRat[1]=surfLib[surf][3]/surfLib[surf][1];	// Omegaphi / Omegar
	OmRat[2]=surfLib[surf][3]/surfLib[surf][2];	// Omegaphi / Omegaz

	/* interpolation */
	/*
	struct gslInterpFunctionParams par;
	par.nn=nn; par.x=&EInterp[0]; par.y=&diffInterp[0];

	gsl_function f;
	f.function = &gslInterpFunction;
	f.params = &par;


	double Enot=0;
	if (surf==0 || surf==nsurfs-1) Enot=-E;
	else {
		gsl_set_error_handler_off();
		const gsl_root_fsolver_type *T;
		gsl_root_fsolver *s;

		int status=0,iter = 0;
		double low=-surfLib[surf-1][0],high=-surfLib[surf+1][0];
		T = gsl_root_fsolver_brent;
		s = gsl_root_fsolver_alloc (T);
		if (gsl_root_fsolver_set (s, &f, low, high)!= GSL_SUCCESS) Enot=-E;
		else {
			while (status == GSL_CONTINUE || high-low>1e-4) {
				iter++;
				status = gsl_root_fsolver_iterate (s);
				Enot = gsl_root_fsolver_root (s);
				low  = gsl_root_fsolver_x_lower (s);
				high = gsl_root_fsolver_x_upper (s);
				status = gsl_root_test_interval (low, high,
												 0, 0.0001);
			}

			gsl_root_fsolver_free (s);
		}

	*/


	if(E<evpot(ar[0],ar[0]) || E>0 || isnan(E)==1) {
		//printf("J=%f %f %f, Efit=%f Exv=%f\n",Jr,Lz,Jz,E,Etrue);
		//fprintf(fdf2,"%e %e %e %e\n",E,Etrue,H_iso(Jr,Lz,Jz),diff);
		return (Etrue>0 ? 0 : Etrue);
	}
	else {//fprintf(fdf2,"%e %e %e %e\n",E,Etrue,H_iso(Jr,Lz,Jz),diff);
		return E;}
}
/*
 *   homogeneous function of degree 1 in (Jr,Lz,Jz)
 */
double hj(double Jr, double Lz, double Jz, double R, double J0){
		double kappa,nu,Omegac;
		double Rc,Lc,Lc_old=0.;
		count++;

		ELz_grid grid;
		Lc=Jr+sqrt(2.)*Lz+sqrt(2.)*Jz;

		//printf("\n>> R=%f, %e %e %e Omegas %e %e\n",R,Jr,Lz,Jz,uvorb.Omegau,uvorb.Omegav);
		int nn=0;
		while(fabs(Lc-Lc_old)>1e-3){
			Lc_old=Lc;

			Rc=grid.GetRc_ini(Lc*Lc,R);
			getfreqs_ini(Rc,&kappa,&nu,&Omegac);
			Lc=Jr+Omegac/kappa*Lz+nu/kappa*Jz;
			if (isnan(Lc)==1) printf("Lc is NANANANANANAN\n");
			nn++; if (nn>50) printf("Difficulty in computing h(J)... Lc=%f Lc_old=%f\n",Lc,Lc_old);
		}
		//fprintf(fdf,"%lf %lf %lf %lf %lf\n",Lz,pow(Lc,-2),R,nu/kappa,Omegac/kappa);


		if (kappa==0.) return -1;
		//return (Jr+J0)+.4*Omegac/kappa*(Lz+J0)+1.2*nu/kappa*(Jz+J0);;
		//	spherical
		//return Lc=Jr+J0+Omegac/kappa*(Lz+J0)+nu/kappa*(Jz+J0);
		//return Jr+fabs(Lz)+Jz;

		//	azimuthal bias deltaphi=.4,deltaz=1.2
		//  radial bias deltaphi=1.4,deltaz=1.5
		double deltar,deltaphi=.4,deltaz=1.2;
		//deltar=1.-Omegac/kappa*(deltaphi-Omegac/kappa)-nu/kappa*(deltaz-nu/kappa);
		deltar=1.-pow(Omegac/kappa,1)*(deltaphi-1.)-pow(nu/kappa,1)*(deltaz-1.);

		Lc=deltar*Jr+deltaphi*Lz+deltaz*Jz;
		//return Lc;

		//  modification for low ellipticities
		//double deltazero=(pow(Omegac/kappa,2)+pow(nu/kappa,2)-deltaz*nu/kappa+1.)/(Omegac/kappa+1.);
		//double deltazero=1.-pow(nu,2)/(pow(Omegac,2)+pow(kappa,2))*(deltaz-1.);
		double deltazero=1-(deltaz-1)/(1+kappa/Omegac);
		double psi=tanh(Lc/sqrt(GM*a));
		deltaphi=(1.-psi)*deltazero+psi*deltaphi;
		deltar=(1.-psi)*deltazero+psi*deltar;

		// use for spherical
		//deltaz=deltaphi;

		Lc=deltar*(Jr+J0)+deltaphi*Omegac/kappa*(Lz+J0)+deltaz*nu/kappa*(Jz+J0);
		//*Om=Omegac/kappa;

		return Lc;
}
/*
 *  f(J)=f[H(J)] \propto h(J)^-2, where h(J)=Jr+Omega_phi/Omega_r*Jphi+Omega_z/Omega_r*z
 *  Isothermal model
 */
double hj_isoth(double Jr, double Lz, double Jz, double R, double Etrue){


	//double mass=3.8e0, J0=1e-1, JM=-evpot(ar[0],ar[0]);

	/*
	 * 18/05/14
	 * LP edit: added surfaces of const E computation.
	 */
	//double j,HJ = EtJ->E(Jr,Jz,Lz,&j,surf_LIB[0]);

	//double OmegaRatios[3]={0.,0.,0.},HJ = quadE(Jr,Lz,Jz,surf_LIB,OmegaRatios,nsurf,Etrue);

	//double j; printf("%f %f %f, E= %f %f\n",Jr,Lz,Jz,HJ,EtJ->E(Jr,Jz,fabs(Lz),&j,surf_LIB[0]));

	//return exp(evpot(ar[0],ar[0])-HJ)/mass;
	//return MAX(0.,pow(1+fabs(evpot(ar[0],ar[0])-HJ),-2)/mass-pow(JM,-2)/mass);
	//if (HJ==0.)return 0; else return MAX(0,pow((J0+Jr)+OmegaRatios[1]*(J0+Lz)+OmegaRatios[0]*(J0+Jz),-2)/mass); //MAX(0,exp(evpot(ar[0],ar[0])-(Jr+OmegaRatios[1]*Lz+OmegaRatios[0]*Jz)))/mass;


	/*
	 * standard version
	 */
	double mass=4e2, J0=1e-1, JM=11.15;
	//double mass=3e4, J0=1e-2, JM=21.15;
	/*
	double J=hj(Jr,Lz,Jz,R,0.);
	//if (J<J0) return pow(J+J0,-1)/mass;//printf("%e %e\n",MAX(0.,pow(J0,-2)/mass-pow(3*JM+J0,-2)/mass),pow(TPI,-3/2)*(exp(H)-1.));
	if (isnan(J)==1 || J<0) return 0;
	return MAX(0.,pow(J+J0,-2)/mass-pow(3*JM,-2)/mass);
	*/

	// hernquist
	double J=Jr+.5*(Lz+4.*Jz); mass=1.4e1; J0=1;
	//fprintf(fdf,"%e %e %e %e\n",MAX(0.,pow(J+J0,-5.)/mass-pow(J0+3*JM,-5.)/mass),Jr,(fabs(Lz)+Jz),Etrue);

	//return MAX(0.,pow(J+J0,-6.)/mass-pow(J0+3*JM,-6.)/mass);
	/*
	 * 10/10/14: adding g(J)!=h(J) in the most general form
	 *           f ~ [J0+hJ]^beta / [ hJ^beta [J0+gJ]^alpha ]
	 *
	 */
	// std hernq
	//return MAX(0.,pow(J,-5./3.)*pow(J0+J,-5+5./3.)/mass);

	// mod henq
	double hJ=Jr+.55*(Lz+Jz),gJ=Jr+1.*(Lz+Jz); mass=1e2;
	return MAX(0.,pow(J0+hJ,5./3.)*pow(hJ,-5./3.)*pow(J0+gJ,-5.)/mass);

	// NFW
	//mass=8e1;
	//return MAX(0.,pow(J,-5/3)*pow(J0+J,-3+5/3)/mass);
}
static double eqApoPeri(double *par, double r){
	double E=par[0],L=par[1];
	return .5*pow(L/r,2.0)+(PhiR(r)-E);
}
static double eqRmax(double *E,double r){
	return PhiR(r)-(*E);
}
double glE,glL;
static double JrInt(double r){
	return sqrt(2.*glE-2.*PhiR(r)-pow(glL/r,2));
}
double df(double *x,double *v){
	double R=x[0],Lz=R*v[1],Phigl=x[2];
	double j;

	/* Hamiltonian: here is needed to compute Delta */
	double H=.5*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2))+Phigl;
	if (H>=0) return 0.;

	/* get the actions */
	double p[2]={v[0],v[2]},Jr,Jz;
	uv_orb uvorb(Deltafn(H),Lz,Phigl,x,p);
	//int nt=int_acts(Lz,uvorb.E,uvorb.Er,&Jr,&Jz);
	//if(nt==0){
		Jr=uvorb.Ju(); Jz=uvorb.Jv();// Jz=JzI;
	//}
	//double z=x[1],phi=0.,L=sqrt( pow(phi*v[2]-z*v[1],2)+pow(z*v[0]-R*v[2],2)+pow(R*v[1]-phi*v[0],2) );
	//printf("Jz+|Lz|=%f L=%f\n",Jz+fabs(Lz),L);

	/*
	double apoperi[2]={0.,0.},par[2],rmax=zbrent(&H,&eqRmax,ar[0],1e10,eqRmax(&H,ar[0]),eqRmax(&H,1e10),1e-5,50);
	par[0]=H; par[1]=L;
	if (eqApoPeri(par,ar[0])*eqApoPeri(par,rmax)<0) apoperi[0]=zbrent(par,&eqApoPeri,ar[0],rmax,eqApoPeri(par,ar[0]),eqApoPeri(par,rmax),1e-5,50);
	else if (eqApoPeri(par,ar[0])*eqApoPeri(par,rmax/2)<0) apoperi[0]=zbrent(par,&eqApoPeri,ar[0],rmax/2,eqApoPeri(par,ar[0]),eqApoPeri(par,rmax/2),1e-5,50);
	else if (eqApoPeri(par,rmax/2)*eqApoPeri(par,rmax)<0) apoperi[0]=zbrent(par,&eqApoPeri,rmax/2,rmax,eqApoPeri(par,rmax/2),eqApoPeri(par,rmax),1e-5,50);

	if (eqApoPeri(par,apoperi[0]*.999)*eqApoPeri(par,ar[0])<0) apoperi[1]=zbrent(par,&eqApoPeri,ar[0],apoperi[0]*.999,eqApoPeri(par,ar[0]),eqApoPeri(par,apoperi[0]*.999),1e-5,50);
	if (eqApoPeri(par,apoperi[0]*1.001)*eqApoPeri(par,rmax)<0) apoperi[1]=zbrent(par,&eqApoPeri,apoperi[0]*1.001,rmax,eqApoPeri(par,apoperi[0]*1.001),eqApoPeri(par,rmax),1e-5,50);
	glE=H; glL=L;
	double Jr2=1./acos(-1.)*oct_int(&JrInt,( apoperi[0]<apoperi[1] ? apoperi[0] : apoperi[1] ),( apoperi[0]<apoperi[1] ? apoperi[1] : apoperi[0] ),12);
	if (isnan(Jr2)!=1) fprintf(fdf2,"%f %f %f %f %f %f\n",Jz+fabs(Lz),L, Jr,Jr2,Jr_iso(H,L),Jr_iso(H,Jz+fabs(Lz)));
	*/

	//uvorb.GetFreqs();
	//double Omegar=uvorb.Omegau,Omegaphi=uvorb.Omegaphi;

	/* compute \bar{j}, here called j */
	//j = (Jr + fabs(Lz) + Jz)/3.;
	j = (Jr + 1.*(fabs(Lz) + Jz))/3.;
	if(j<0) printf("error: j<0 %g ",j);

#if defined(HERNQUIST) || defined(NFW)			// use the same procedure for Hernquist and NFW models
	return hj_isoth(Jr,Lz,Jz,R,H);
	return expdf_H(j,Jr,Lz,Jz,x,v,5e-2,-1,5e-2);
			//+ .5*expdf_H(j,Jr,Lz,Jz,x,v,1.e-5,1e-5,-1.);

	/* Pontzen & Governato */
	//return 	expdf_H(j,Jr,Lz,Jz,x,v,2.86e-4,1.69e-3,1.59e-3);
#elif defined ISOTHERMAL

	return hj_isoth(Jr,Lz,Jz,R,H);

#elif defined ISOCHRONE
	//	Jr=Jr_iso(H,L); Jz=L-fabs(Lz);
	return hj_isoth(Jr,Lz,Jz,R,H);
	return df_I(Jr,Lz,Jz,R,H);
	//return expdf_I(j,Jr,Lz,Jz,R,0,0,0,H);

	/* Pontzen & Governato */
	//return 	expdf(j,Jr,Lz,Jz,2.86e-4,1.69e-3,1.59e-3);
#endif
}

double rho_iso(double r){
	double t=sqrt(r*r+a*a);
	return GM*(3*(a+t)*t*t-r*r*(a+3*t))/pow((a+t)*t,3)*.5/TPI;
}
double phi_iso(double r,double b){
	return -GM/(b+sqrt(b*b+r*r));
}
