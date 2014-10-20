/*
 * surf_lib.cpp
 *
 *  Created on: Mar 12, 2013
 *      Author: Fermani
 */

#include "surf_lib.h"
#include "Types.h"
#include "Pi.h"
#include "../oct_int_exp.h"
#include "tools.h"
#include "eqs_of_motion.h"
#include "Iso_check.h"
#include "LM_fit.h"
#include "../press.h"
#include "computeJtheta.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#define ITMAX 100
#define EPS 3.0e-8

using namespace std;

using std::cout;
using std::cerr;

eveq_p gl_par; // needed for Jr_integrand
int info; // to require extra info in debugging
int verbose=1;

surf_lib::surf_lib(eveq_p par,string file,int length){
	length_=length;
	a_=par;
	output_=file;
	list_ = new double* [length];
	for(int i=0;i!=length;i++) list_[i] = new double[10];
}

surf_lib::~surf_lib(){
	for(int i=0;i!=length_;i++) delete[] list_[i];
	delete[] list_;

	delete &length_;
	delete &a_;
	delete &output_;
}

/*
 * 12/05/14
 * LP edit: overloaded member build(). The function will initialize
 * 			the energy limits.
 */
void surf_lib::build(){
	this->build(1,.5);
}

double RcEfn(double *Ep,double R){//what vanishes at Rc for given E
	return .5*R*gl_par.p->dphidR(R,0)+gl_par.p->phi(R,0)-(*Ep);
}
double RcE(double E,double Rc){
	double R0=.8*Rc,f0=RcEfn(&E,R0);
	while(f0>0){
		R0*=.8; f0=RcEfn(&E,R0);
	}
	double R1=1.2*Rc, f1=RcEfn(&E,R1);
	while(f1<0){
		R1*=1.2;  f1=RcEfn(&E,R1);
	}
	return zbrent(&E,&RcEfn,R0,R1,f0,f1,1.e-5,25);
}

void surf_lib::build(double minE, double maxE){
	vec2 find_apo_peri(eveq_p *par);
	double Jr_integrand(double r);
	double find_R(eveq_p *par,double eps,int* status);
	double find_Rg(eveq_p *par);

	/*
	 * 12/05/14
	 * LP edit: energy limits as parameters for build
	 */
	if (maxE>0.) maxE = a_.p->phi(0.05,0.05);
	if (minE>0.) minE = a_.p->phi(50,50);
	if (maxE > minE){
		cerr << "Maximum (abs. value) requested for the energy is greater than "
			 <<	"the minimum (abs. value)...stopping!";
		cout << "maxE=" << a_.p->phi(0.,0.) << ", minE=" << a_.p->phi(50.,50.) << endl ;
		exit(1);
	}

	cout << maxE << "\n";
	double eps_e=1.0e-1, eps=1.0e-6,E;
	double Jr,Jz,Jphi;
	double Rg,maxJphi,maxJz;
	vec2 apo_peri; eveq_p par=a_;
	int i,status=0,k,ite,jl; // if SUCCESS, status=0

	//test of surf_build;
	ofstream to;
	string output("outputB.dat");
	const char *o=(const char*)output.c_str();
	Iso_check *tiso;
	double coefficients[2]={par.p->M_coeff(), par.p->b_coeff()};
	tiso=new Iso_check(2,&coefficients[0]);
	int test=0;
	if(test==0) to.open(o);
	double temp_add,perturb_E=0.;

	double **Jlist= new double* [21];
	for(i=0;i<21;i++) Jlist[i]= new double[3];

	LM_fit *fit;
	fit = new LM_fit(par.energy,21,Jlist);

	to << scientific;

	i=0;
	while(i<length_){
		if(verbose) cout << endl << " ==== starting: " << i << "/" << length_ << endl;
		test=i;
		//E=maxE*double(i)/double(length_-1)-eps_e-perturb_E;
		//E=(maxE-minE)*double(i)/double(length_)+minE-perturb_E;
		E=maxE-sinh(i*asinh(maxE)/(double)(length_-1)) -perturb_E; if (i==length_-1) E+=minE;

		/*
		 * LP EDIT 30/05/14: removed this if
		 *
		if(fabs(E)>fabs(maxE)){
			temp_add=fabs((maxE*double(i)/double(length_-1)-eps_e)-(maxE*double(i-1)/double(length_-1)-eps_e));
			while(fabs(E)>fabs(maxE)){
				temp_add/=1.5;
				E=maxE*double(i-1)/double(length_-1)-eps_e-temp_add;
			}
			if(verbose) cout << "Energy adapted to " << E << "\n";
		}
		*/
		par.energy=E;
		fit->reset_energy(E);
		jl=0;

		//SURFACE CORNERS:
		// Jr 0 0 point
		info=0;
		Jz=0.; Jphi=0.;
		par.Lz=Jphi;
		par.L=Jz+Jphi;
		apo_peri=find_apo_peri(&par);
		gl_par=par; //initialization needed by Jr_integrand
		Jr=1.0/Pi*oct_int(&Jr_integrand,apo_peri[0],apo_peri[1],lmax);
		Jlist[jl][0]=Jr; Jlist[jl][1]=Jz; Jlist[jl][2]=Jphi; jl++;
		if(test==0 || info==1){
			if(verbose) cout << "info " << info << "\n";
			to << Jr << " " << Jz << " " << Jphi << " " << tiso->Jr(E,Jz,Jphi) << " " << Jz << " " << Jphi << " " "\n";
			if(verbose) cout << Jr << " " << Jz << " " << Jphi << " " << tiso->Jr(E,Jz,Jphi) << " " << Jz << " " << Jphi << " " "\n";
		}

		// 0 Jz 0
		info=0;
		Jr=0.; Jphi=0.;
		Jz=find_R(&par,eps,&status); //modify input, make output be Jz and get it to modify status to 1 if failure. here there is orb computation, need to insert status cnd
		if(status!=0){
			ite=0;
			while(ite<3 && status!=0){
				if(verbose) cerr << "(0,Jz,0) point: computation failed. \n Try to adjust integration tolerance";
				if(verbose) cout << "[Jz,status]=[" << Jz << "," << status << "]\n";
				eps=eps*10.0;
				Jz=find_R(&par,eps,&status);
				ite++;
			}
			if(status!=0){
				if(verbose) cerr << "(0,Jz,0) point: computation failed. \n Try to adjust energy lower bound (all points reseted)";
				eps_e=eps_e*2.0;

			}
		}
		maxJz=Jz;
		Jlist[jl][0]=Jr; Jlist[jl][1]=Jz; Jlist[jl][2]=Jphi; jl++;
		if(test==0 || info==1){
			if(verbose) cout << "info " << info << "\n";
			to << Jr << " " << Jz << " " << Jphi << " " << Jr << " " << tiso->L(E,Jr)-Jphi << " " << Jphi << " " "\n";
			if(verbose) cout << Jr << " " << Jz << " " << Jphi << " " << Jr << " " << tiso->L(E,Jr)-Jphi << " " << Jphi << " " "\n";
		}

		if(status==0){
			eps=1.0e-6;
			// 0 0 Jphi
			info=0;
			Jr=0.; Jz=0.;
			//Rg=find_Rg(&par);
			Rg=RcE(E,.1);
			Jphi=pow(Rg,1.5)*sqrt(par.p->dphidR(Rg,0.0));
			//cout << "RG: " << Rg   << " RCE: " << RcE(E,1) << endl;
			//cout << "RG: " << Jphi << " RCE: " << pow(RcE(E,1),1.5)*sqrt(par.p->dphidR(RcE(E,1),0.0)) << endl;
			maxJphi=Jphi;
			Jlist[jl][0]=Jr; Jlist[jl][1]=Jz; Jlist[jl][2]=Jphi; jl++;
			if(test==0 || info==1){
				if(verbose) cout << "info " << info << "\n";
				to << Jr << " " << Jz << " " << Jphi << " " << Jr << " " << Jz << " " << tiso->L(E,Jr)-Jz << " " "\n";
				if(verbose) cout << Jr << " " << Jz << " " << Jphi << " " << Jr << " " << Jz << " " << tiso->L(E,Jr)-Jz << " " "\n";
			}

			// Shell orbits (0,x,J_phi)
			for(k=0;k<6;k++){
				info=0;
				Jr=0.;
				eps=1.0e-6;
				Jphi=maxJphi*(0.125*double(k+1));
				par.Lz=Jphi;
				//Jz=maxJz-Jphi;

				Jz=find_R(&par,eps,&status);
				if(status!=0){
					ite=0;
					while(ite<3 && status!=0){
						if(verbose) cerr << "Shell orbits computation failed. \n Try to adjust integration tolerance";
						eps=eps*10.0;
						Jz=find_R(&par,eps,&status);
						ite++;
					}
					if(status!=0){ cerr << "Shell orbits computation failed repeatdly. Quit \n"; exit(1);}
				}
				Jlist[jl][0]=Jr; Jlist[jl][1]=Jz; Jlist[jl][2]=Jphi; jl++;
				if(test==0 || info==1){
					if(verbose) cout << "info " << info << "\n";
					to << Jr << " " << Jz << " " << Jphi << " " << Jr << " " << tiso->L(E,Jr)-Jphi << " " << Jphi << " " "\n";
					if(verbose) cout << Jr << " " << Jz << " " << Jphi << " " << Jr << " " << tiso->L(E,Jr)-Jphi << " " << Jphi << " " "\n";
				}
			}

			// Plane orbits (x,0,J_phi)
			for(k=0;k<6;k++){
				info=0;
				Jz=0.;
				eps=1.0e-6;
				Jphi=maxJphi*(0.15*double(k+1));
				par.Lz=Jphi;
				par.L=Jz+Jphi;
				apo_peri=find_apo_peri(&par);
				gl_par=par; //initialization needed by Jr_integrand
				Jr=1.0/Pi*oct_int(&Jr_integrand,apo_peri[0],apo_peri[1],lmax);

				Jlist[jl][0]=Jr; Jlist[jl][1]=Jz; Jlist[jl][2]=Jphi; jl++;
				if(test==0 || info==1){
					if(verbose) cout << "info " << info << "\n";
					to << Jr << " " << Jz << " " << Jphi << " " << tiso->Jr(E,Jz,Jphi) << " " << Jz << " " << Jphi << " " "\n";
					if(verbose) cout << Jr << " " << Jz << " " << Jphi << " " << tiso->Jr(E,Jz,Jphi) << " " << Jz << " " << Jphi << " " "\n";
				}
			}

			// Box orbits (x,Jz,0)
			for(k=0;k<6;k++){
				info=0;
				Jphi=0.;
				eps=1.0e-6;
				Jz=maxJphi*(0.15*double(k+1));
				//Jz=maxJz*(0.15*double(k+1));
				par.Lz=Jphi;
				par.L=Jz+Jphi;
				apo_peri=find_apo_peri(&par);
				gl_par=par; //initialization needed by Jr_integrand
				Jr=1.0/Pi*oct_int(&Jr_integrand,apo_peri[0],apo_peri[1],lmax);

				Jlist[jl][0]=Jr; Jlist[jl][1]=Jz; Jlist[jl][2]=Jphi; jl++;
				if(test==0 || info==1){
					if(verbose) cout << "info " << info << " maxJz= " << maxJz << " max Jphi=" << maxJphi << "\n";
					to << Jr << " " << Jz << " " << Jphi << " " << tiso->Jr(E,Jz,Jphi) << " " << Jz << " " << Jphi << " " "\n";
					if(verbose) cout << Jr << " " << Jz << " " << Jphi << " " << tiso->Jr(E,Jz,Jphi) << " " << Jz << " " << Jphi << " " "\n";
				}
			}
			fit->fit_qsurface_2list(100);
			if(fit->status()==0){
				for(k=0;k<10;k++) list_[i][k]=0.;
				list_[i][0]=E;
				for(k=1;k<=fit->ncoeff();k++){ list_[i][k]=fit->surf_coeff(k-1); cout << fit->surf_coeff(k-1) << " ";}
				if(verbose) cout << "\n";
				perturb_E=0.0;
				i++; if(!verbose) cout << i << "/" << length_ << ", ";
			}
			else{
				perturb_E+=0.0005*E; // if unable to fit surface, shift the energy by 0.5% recompute.
				if(verbose) cerr << "Unable to fit surface, try to shift energy by 0.05% and recompute.\n";
			}
		}
		if(test==0){ if(verbose) cout << "Energy of surf. edges compared: " << E << "\n"; to.close();}
		else if(verbose) cout << "E surf. edges computed: " << E << "\n";

		//exit(1);
	}
	if(verbose) cout << "MaxE MaxE_obtained: " << maxE << " " << E << "\n";
}

void surf_lib::save(double **surf_LIB){
	const char *o=(const char*)output_.c_str();
	ofstream to;
	to.open(o);
	for(int i=0;i<length_;i++){
		for(int k=0;k<9;k++){
			surf_LIB[i][k]=list_[i][k];
			to << list_[i][k] << " ";
		}
		surf_LIB[i][9]=list_[i][9];
		to << list_[i][9] << "\n";
	}
	to.close();
}


int surf_lib::read(double **surf_LIB){
	const char *in=(const char*)output_.c_str();
	FILE *trojab = fopen(in,"r");
	char buffer2[8192];
	char **numen = NULL;
	int nzb = 0, k;
	while (fgets(buffer2, 8191, trojab)) {
		int nnb = 0;
		for(k=0;k<10;k++)
		{
			surf_LIB[nzb][k] = strtod(&buffer2[nnb], numen);
			while (buffer2[nnb] != ' ') {
				nnb++;
			}
			while (buffer2[nnb] == ' ') {
				nnb++;
			}
		}
		nzb++;
	}
	fclose(trojab);
	return nzb;
}

vec2 find_apo_peri(eveq_p *par){
	double zbrent(double eq_sm,double at_this_z,eveq_p *par,double (*func)(double,double,eveq_p*), double x1, double x2, double tol);
	double eq_apo_peri(double r,double z0,eveq_p *par);
	double compute_rmax(eveq_p *par);
	vec2 res;

	double temp,RMAX=compute_rmax(par);
	//cout << "Rmax " << RMAX << "\n";

	res[0]=zbrent(0.,0.,par,&eq_apo_peri,rinf,RMAX,tollRF);
	// Specific for eq_apo_peri: (the function within the roots is negative) //
	if(eq_apo_peri(res[0]*.999,0.,par)<eq_apo_peri(res[0],0.,par)) res[1]=zbrent(0.,0.,par,&eq_apo_peri,rinf,res[0]*.998,tollRF);
	else res[1]=zbrent(0.,0.,par,&eq_apo_peri,res[0]*1.002,RMAX,tollRF);
	//cout << "apo_peri [1]" << res[0] << " " << res[1] << "\n";
	if(fabs(res[0]-res[1])<=min(res[0],res[1])*0.004){
		if(par->L!=0){
			if(verbose) cerr << "You are considering orbit going through the center, but L!=0..mmm\n";
			if(verbose) cout << "apo_peri [1]" << res[0] << " " << res[1] << "\n";
			info=1;
		}
		if(res[0]>res[1]) res[1]=rinf;
		else res[0]=rinf;
	}

	if(res[0]>res[1]){ temp=res[0]; res[0]=res[1]; res[1]=temp;}
	//cout << "apo_peri [2]" << res[0] << " " << res[1] << "\n";
	return res;
}


double eq_apo_peri(double r,double z0,eveq_p *par){
	z0=0.;
	if(par->L<EPS && r<EPS) return 0.;
	else return pow(par->L/r,2.0)+2.0*(par->p->phi(r,z0)-par->energy);
}

double compute_rmax(eveq_p *par){
	double compute_Rmax_eff(double at_z,eveq_p *par);
	return compute_Rmax_eff(0.0,par);
}

double compute_Rmax_eff(double at_z,eveq_p *par){
	double zbrent(double eq_sm,double at_this_z,eveq_p *par,double (*func)(double,double,eveq_p*), double x1, double x2, double tol);
	double eq_for_Rmax(double R,double z,eveq_p *par);

	if(par->energy>=par->p->phi(rsup,at_z)){
		cerr << "Energy too close to zero. sup(r) exceeds maximum value: please change rsup in Structure.h \n";
		cout << "Orbits will go further than: " << rsup << " kpc.\n";
		cout << par->energy << " " << par->p->phi(rsup,at_z) << "\n";
		exit(1);
	}
	return zbrent(par->energy,at_z,par,&eq_for_Rmax,rinf,rsup,tollRF);
}

double eq_for_Rmax(double R,double z,eveq_p *par){
	return par->p->phi(R,z);//+0.5*pow(par->Lz/R,2.0);
}

double Jr_integrand(double r){
	if(fabs(2.0*gl_par.energy-2.0*gl_par.p->phi(r,0.)-pow(gl_par.L/r,2.0))<1e-5) return 0.0;
	/*
	 * 05/06/14
	 * LP edit: I get nans sometimes, so I set return 0 in those cases...
	 */
	else return (isnan(sqrt(2.0*gl_par.energy-2.0*gl_par.p->phi(r,0.)-pow(gl_par.L/r,2.0)))==1 ?
				 0:
				 sqrt(2.0*gl_par.energy-2.0*gl_par.p->phi(r,0.)-pow(gl_par.L/r,2.0)));
}

double find_Rg(eveq_p *par){
	double compute_Rmax_eff(double at_z,eveq_p *par);
	double zbrent(double eq_sm,double at_this_z,eveq_p *par,double (*func)(double,double,eveq_p*), double x1, double x2, double tol);
	double guiding_centre_eq(double R,double z,eveq_p *par);

	//cout << "Par passage check (at find_Rg): " << params->energy << " " << params->p_ax_pot(0.0,0.0) << "\n";
	double supRg=compute_Rmax_eff(0.0,par);
	return zbrent(par->energy,0.,par,&guiding_centre_eq,0.0,supRg,toll);
}


double guiding_centre_eq(double R,double z,eveq_p *par){ return 0.5*R*par->p->dphidR(R,z)+par->p->phi(R,z);}

double zbrent(double eq_sm,double at_this_z,eveq_p *par,double (*func)(double,double,eveq_p*), double x1, double x2, double tol)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a,at_this_z,par)-eq_sm,fb=(*func)(b,at_this_z,par)-eq_sm,fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){}
		//cerr << "Root must be bracketed in zbrent \n";
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += sign(tol1,xm);
		fb=(*func)(b,at_this_z,par)-eq_sm;
	}
	cerr << "Maximum number of iterations exceeded in zbrent \n";
	return 0.0;
}
void checkImagPz(eveq_p * par, double launch_R){
	if(2.0*(par->energy-par->p->phi(launch_R,0.0))-pow(par->Lz/launch_R,2.0)<0){
		/*
		cerr << "in find_R, pz is imaginary. (R,z,pR,pz)=(" << launch_R << ",0,0,sqrt(" << 2.0*(par->energy-par->p->phi(launch_R,0.0)-0.5*pow(par->Lz/launch_R,2.0)) << "))\n";
		cout << "E " <<  par->energy << "\n";
		cout << "pot " << par->p->phi(launch_R,0.0) << "\n";
		cout << "Lz term " << pow(par->Lz/launch_R,2.0) << "\n";
		*/
		throw 1;
	//exit(1);
	}
	else throw 0;
}
double find_R(eveq_p *par,double eps,int* status){
	double compute_Rmax_eff(double at_this_z,eveq_p *par);
	double sup_R, inf_R;
	double compute_R_crossing_meridial_plane(double init_R,double p_z,eveq_p *params,double *Jz, int *status,double eps);
	double launch_R,land_R, pz, Jz;
	int iteration=0;

	double rg=RcE(par->energy,.1);
		par->L=pow(rg,1.5)*sqrt(par->p->dphidR(rg,0.0));
		//par->Lz=pow(rg,1.5)*sqrt(par->p->dphidR(rg,0.0))/5;
		cout << "E,L,Lz: " << par->energy << " " << par->L << " " << par->Lz << endl;
		Jz= getJtheta(par,rg);	//2* for isochrone??


		//exit(1);
		return Jz;


	sup_R=compute_Rmax_eff(0.0,par);
	/*
	 * 15/05/14
	 * LP edit: changed inf_R=0 to inf_R=0.01, since my potentials stop before r=0
	 */
	inf_R=0.005;

	launch_R=0.5*(sup_R-inf_R)+inf_R;
	land_R=launch_R+10.0*toll;
	//cout <<  sup_R << ' ' << ' ' << inf_R << endl;
	while(fabs(launch_R-land_R)>toll && iteration<max_ite){
		launch_R=0.5*(sup_R-inf_R)+inf_R;

		/*
		 * 03/06/14
		 * LP Edit: exception handling added
		 */
		int b=0;
		while(b<1){
			try {
				checkImagPz(par,launch_R);
			}
			catch (int e){
				if (e==0) b=1;
				else launch_R *= 1.025;
			}
		}
		pz=sqrt(2.0*(par->energy-par->p->phi(launch_R,0.0)-0.5*pow(par->Lz/launch_R,2.0)));

		land_R=compute_R_crossing_meridial_plane(launch_R,pz,par,&Jz,status,eps);
		if(verbose) if(*status!=0) cerr << "status for land_R is " << *status << ". [E,pz,launch_R,land_R]=[" << par->energy << "," << pz << "," << launch_R << "," << land_R << "]. Ite n." << iteration << "\n";
		//cout << " -- " << iteration << ' ' << land_R << endl;
		if(fabs(launch_R-land_R)<toll) iteration=ID+1;
		else{
			if(land_R>launch_R) inf_R=launch_R;
			else sup_R=launch_R;
			iteration++;
		}
	}
	/*
	 * Try integrating directly: Jz = 1/Pi int_0^zmax sqrt(2*(E-Phi(Rc,z))) dz
	 */

	if(fabs(launch_R-land_R)>toll && iteration==max_ite){
		if(verbose) cerr << "When computing z-orbits, iteration for pz estimate failed. \n. Energy=" << par->energy << ", [R_launch,R_land]=[" << launch_R << "," << land_R << "] kpc.\n";
		//*status=1;
	}

	double computeJz(eveq_p * par, double Rc);
	//cout << "Jz COMPARISON: " << Jz << " " << computeJz(par,land_R) << endl;


	//exit(1);
	return Jz;
}

double compute_R_crossing_meridial_plane(double init_R,double p_z,eveq_p *par,double *Jz, int *status,double eps){
	double res=0.0, ti=0.0;
	vec2 pA,pB;
	bool condition=0;
	*status=0;

	// GSL RK8
	double y[4] = {init_R, 0.0, 0.0, +p_z };
	if (isnan(p_z)==1) {
		cout << "Pz=NAN!" << endl;
		exit(1);
	}

	gsl_odeiv2_system sys = {eqs_of_motion, NULL, 4, par}; // GSL integrator tested and functioning
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, eps, eps, 0.0);


	double t=0.0;
	double grad=1.0, intercept=0.0;
	double mR=rsup;
	double Mz=0.0;
	pA[0]=init_R; pA[1]=0.0; pB=pA;

//	to.open(out_p);
	int stat,i=0;

	vec2 orbitA, orbitB;
	orbitB[0]=0.0; orbitB[1]=p_z;
	*Jz=0.0;
	while(condition==0 && i<(int)(ID/step)){
		orbitA=orbitB;
		pA=pB;
		//cout << t << " " << y[1] << " " << y[2] << " " << y[3] << " pz=" << p_z << endl;
		stat = gsl_odeiv2_driver_apply (d, &ti, t, y);
		if (stat != GSL_SUCCESS){
			printf ("error, return value=%d\n", stat);
			if(verbose) cout << "[E,R,pz]=[" << par->energy << " , " << init_R << " , " << p_z << "]\n";
			condition=1;
			*status=1;
			exit(1);
		}
		orbitB[0]=y[1]; orbitB[1]=y[3];

		//cout << "t= " << t << " R= " << y[0] << " z= " << y[1] << " pR= " << y[2] << " pz= " << y[3] << endl;

		pB[0]=y[0]; pB[1]=y[1];

		if(y[1]<0 && i>0){
			condition=1;
			grad=(pA[1]-pB[1])/(pA[0]-pB[0]); // linear interpolation to find R at z=0;
			intercept=pA[1]-grad*pA[0];
			res=-intercept/grad;

			grad=(orbitA[1]-orbitB[1])/(orbitA[0]-orbitB[0]); // linear interpolation to find pz at z=0;
			intercept=orbitA[1]-grad*orbitA[0];
			orbitB[0]=0.0; orbitB[1]=intercept;
			*Jz+=(orbitA[1]+orbitB[1])*(orbitB[0]-orbitA[0])/2.0;

			grad=(pB[1]-pA[1])/(step); // linear interpolation to find t at z=0;
			intercept=pB[1]-grad*t;
			//Omegaz_freq=Pi/fabs(-intercept/grad);
		}
		else{
			if(fabs(pB[0])<mR) mR=fabs(pB[0]);
			if(pB[1]>Mz) Mz=pB[1];
			*Jz+=(orbitA[1]+orbitB[1])*(orbitB[0]-orbitA[0])/2.0;
			//Omegaphi_freq+=GL_par.Lz/pow(pB[0],2.0);
			i++;
		}
		t=t+step/100.0;
	}

	if(condition==0 && i>=(int)(ID/step)){
		//cerr << "Orbit integrated for more than " << ID << " time steps. Are you sure this is right?! I DON'T! \n";
		*status=1;
	}

	*Jz=1.0/(Pi-acos(mR/sqrt(mR*mR+Mz*Mz)))*(*Jz);
	gsl_odeiv2_driver_free (d);
	return fabs(res);
}

struct JzIntPar {int len; double *x; double *y;};
double JzIntegrand(double z, void* params){

	struct JzIntPar * par = (struct JzIntPar *) params;
	int len = (par->len);
	double *x = (par->x);
	double *y = (par->y);

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline,len);
	gsl_spline_init(spline, x, y, len);

	double pz = gsl_spline_eval(spline, z, acc);

	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

	return pz;
}

double evolveOrbit(eveq_p *par, double *zmax, int *n, double *Rini, double *z,double *Pz, double *R, int size, double steph) {
	double pz,y[4],eps=1e-6;
	int i,nn,k=0; double t,ti;

	for (int j=0; j<size; j++) {Pz[j]=0; z[j]=0; R[j]=0;}
	t=0; ti=0;

	int b=0,c=0;
	while (b<1){
		try{
			pz = sqrt(2.0*(par->energy-par->p->phi((*Rini),0.0)-0.5*pow(par->Lz/(*Rini),2.0)));
			if (isnan(pz)==1) throw 1;
			b=1; c=1;
		}
		catch (int e){
			k++;
			if (e==1) {if (c<1) (*Rini)*=1.025; else (*Rini)*=.975;}
			if (k>100) return 2*(*Rini);
		}
	}

	/*
	y[0]= (*Rini); y[1]= 0.0; y[2]= 0.0; y[3]= +pz;
	gsl_odeiv2_system sys = {eqs_of_motion, NULL, 4, par}; // GSL integrator tested and functioning
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, steph, eps, eps);
	*/

	double L = pow(RcE(par->energy,.1),1.5)*sqrt(par->p->dphidR(RcE(par->energy,.1),0.0));
	y[0]= (*Rini); y[1]= Pi/2; y[2]= 0; y[3]= -(pow(L,2)-pow(par->Lz,2));
	gsl_odeiv2_system sys = {eqs_of_motion_spherical, NULL, 4, par}; // GSL integrator tested and functioning
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, steph, eps, eps);


	*zmax=0; i=0; (*n)=0; nn=0;
	z[0]=0.; Pz[0]=pz; R[0]=(*Rini);

	FILE* orbp = fopen("orbit.dat","w");
	while (i<size && y[1] >= 0.){
		i++;
		gsl_odeiv2_driver_apply (d, &ti, t, y);

		if (y[1]>*zmax) *zmax=y[1];
		Pz[i]=y[3]; z[i]=y[1]; R[i]=y[0];
		ti=t;
		t=t+steph;
		fprintf(orbp,"%f %f %f %f\n",y[0],y[1],y[2],y[3]);
	}
	fclose(orbp);

	// find index of zmax and the meridional crossing Radius
	i=0; while(z[i]<*zmax) {(*n)++; i++;}
	i=0; while(z[i]>=0) {nn++; i++;}

	gsl_odeiv2_driver_free (d);

	/*
	 * linear interpolation to get the meridional crossing radius
	 */
	double m=(z[nn]-z[nn-1])/(R[nn]-R[nn-1]);
	return (z[nn]-m*R[nn])/m;
}

struct eqInversionZParams { eveq_p * par; double R;};
double eqInversionZ(double z, void *params){

	struct eqInversionZParams* p = (eqInversionZParams *) params;
	eveq_p * par = p->par;
	double R = p->R;

	return par->energy-par->p->phi(R,z);
}

double computeZinv(eveq_p * par, double Rc){
	/*
	 * Find approximation of z_inversion with GSL root finder
	 */
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	gsl_function FInvZ;

	struct eqInversionZParams params = {par,Rc};
	FInvZ.function = &eqInversionZ;
	FInvZ.params = &params;

	int status=0,iter = 0;
	double low=0,high=50,zinv=0;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &FInvZ, low, high);

	while (status == GSL_CONTINUE || high-low>1e-3) {
		iter++;
		status = gsl_root_fsolver_iterate (s);
		zinv = gsl_root_fsolver_root (s);
		low  = gsl_root_fsolver_x_lower (s);
		high = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (low, high,
		                                 0, 0.001);
	}

	gsl_root_fsolver_free (s);
	return zinv;
}
double computeJz(eveq_p * par, double Rc){

	double zinv=computeZinv(par,Rc);
	int n=0,k=0; double steph= .01 * zinv / sqrt(fabs(par->energy));
	int size=(int) 1e4; //par->p->phi(.07,.07)/par->energy/steph;

	//cout << size*steph << " " << 1e3*steph << endl; exit(1);

	double Pz[size],z[size],R[size],Rini1=Rc,Rini2=Rc,Rfin=2*Rc;
	double supR = compute_Rmax_eff(0.0,par), infR = 0.001;
	double Rmid, d1=1e2, d2=1e2, zmax, Rini=Rc;

	/*
	 * using Rc as initial condition for the orbit [and returning
	 * 2/(Pi-acos...)] I get the same values as Francesco
	 */
	/*
	while ( d1>toll && d2>toll ){
		k++;
		Rmid = .5 * (supR+infR);

		Rini1=.5*(Rmid+supR);
		Rfin = evolveOrbit(par,&zmax,&n,&Rini1,z,Pz,R,size,steph);
		d1=fabs(fabs(Rini1)-fabs(Rfin));

		Rini2=.5*(Rmid+infR);
		Rfin = evolveOrbit(par,&zmax,&n,&Rini2,z,Pz,R,size,steph);
		d2=fabs(fabs(Rini2)-fabs(Rfin));

		cout << " ## BRACKETING: " << infR << " " << supR << " " << .5*(Rmid+supR) << " " << .5*(Rmid+infR) << " " << d1 << " " << d2 << " LEN: " << n << endl;

		if (d1<d2) {infR=Rmid; Rini=Rini1;}
		else {supR=Rmid; Rini=Rini2;}

		if (k>30) break;
	}
	*/
	Rfin = evolveOrbit(par,&zmax,&n,&Rini,z,Pz,R,size,steph);
	if (n<1) return 0.;

	cout << "FINAL " << Rini << " " << Rfin << " " << fabs(fabs(Rini)-fabs(Rfin)) << " LEN: " << n << "," << k << endl;
	double Pz2[n],z2[n];
	for (int i=0; i<n; i++){
		Pz2[i] = Pz[i+1];
		z2[i]  = z[i+1];
		//cout << z2[i] << " " << Pz2[i] << " " << zmax << endl;
	}

	gsl_set_error_handler_off();
	const gsl_interp_type *interpType;
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	if (n>2) interpType = gsl_interp_cspline;
	else interpType = gsl_interp_linear;

	gsl_interp * interp = gsl_interp_alloc(interpType, n);
	if (gsl_interp_init(interp, z2, Pz2, n) == GSL_EINVAL) {
		cout << "k=" << k << ",  Rini,zm,zi: " << Rini << " " << zmax << " " << zinv << endl;
		for (int i=0;i<n;i++) cout << "z,Pz: " << z2[i] << " " << Pz2[i] << " " << z[i] << " " << Pz[i] << endl;
		exit(1);
	}
	double res = gsl_interp_eval_integ (interp, z2, Pz2, 0, zmax, acc);

	gsl_interp_free (interp);
	gsl_interp_accel_free (acc);

	double mR=R[n];
	cout << "correction: " << mR << " " << acos(mR/sqrt(mR*mR+zmax*zmax))/Pi << " " << mR/sqrt(mR*mR+zmax*zmax) << endl;

	return 2*res/(Pi-acos(mR/sqrt(mR*mR+zmax*zmax)));
}
