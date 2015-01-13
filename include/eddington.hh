/*
 * eddington.hh
 *
 *  Created on: 05/feb/2014
 *      Author: LP
 */
#include <math.h>
#include "press.h"
#include "isodf4.h"
#include "isopot.h"
#include "leg_pot3.h"
#include "progressbar.hh"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_integration.h"

class Edd_DF
{
public:
	Edd_DF(double* R_grid_in, int NN_in);		// Constructor
	~Edd_DF();									// Destructor
	Edd_DF(const Edd_DF &L);       		      	// copy constructor
	Edd_DF & operator=(const Edd_DF &L);	 	// assignment

	double* DF;									// Discretized DF
	double* psi_r;								// Potential grid

private:
	int NN;										// Size of grids
	double Ehere;								// Energy in the integrand
	double *R_grid,*rho_r,*integral,*R_fow;		// radius,radial density, potential and integral grid
	gsl_interp *interpolation;
	gsl_interp_accel * accelerator;
	void create_grid();							// Creates array *integral
	double integrand(double psi);
	double edd_inversion(double E);
	double integrate(double E, int it);
	double dinterp(double *xm,double *ym,
			         int np,double x);
	double trapzd(double a, double b,
				  double *s,int *it, int n);
	double romberg(double E, int it);
	double ddinterp(double *xm, double *ym,
					 int np, double x);
	static double integrand_wrapper(double x,
					 void *params);
};


/*
 *  Constructor: takes the radial grid and its size as inputs
 *  builds the DF array and the psi_r, which gives the energy
 *  where the DF is calculated
 */
Edd_DF::Edd_DF(double* R_grid_in, int NN_in){
	NN=NN_in;
	/* allocations */
	DF=dmatrix(NN); R_grid=dmatrix(NN); R_fow=dmatrix(NN);
	rho_r=dmatrix(NN); psi_r=dmatrix(NN); integral=dmatrix(NN);

	for (int i=0; i<NN; i++) {
		R_grid[i]=R_grid_in[NN-1-i];
		R_fow[i]=R_grid_in[i];
	}

	/* gsl interpolation allocation */
	interpolation = gsl_interp_alloc (gsl_interp_cspline,NN);
    accelerator =  gsl_interp_accel_alloc();

	create_grid();

	/*
	 * eq (4.46a): uses first derivative drho/dpsi
	 */
	 //for (int i=0; i<NN; i++) DF[i]=MAX(0.,edd_inversion(psi_r[i+1]));

	/*
	 * eq (4.46b): uses second derivative d2rho/dpsi2
	 */
	double drhodpsi_zero=rho_r[0]/psi_r[0];
	for (int i=0; i<NN; i++) DF[i]=MAX(0.,1./(sqrt(8.)*M_PI*M_PI)*(integral[i]+drhodpsi_zero/sqrt(psi_r[i])));

	gsl_interp_accel_free (accelerator);
}

/*
 *  Destructor: deletes DF and psi_r
 */
Edd_DF::~Edd_DF(){
	delete [] DF; delete [] psi_r;
	delete [] rho_r; delete [] integral; delete [] R_grid; delete [] R_fow;
}

/*
 *  create_grid: performs the first step of the formula. Integrates the
 *  function integrand(double) and stores the data in the array *integral.
 *  The integral is done with the trapezoidal rule.
 */
void Edd_DF::create_grid(){

#ifdef HERNQUIST
	for (int i=0; i<NN; i++){
		rho_r[i]=rho_Hern(R_grid[i]);
		psi_r[i]=-phi_Hern(R_grid[i]);
	}
#elif defined NFW
	for (int i=0; i<NN; i++){
		rho_r[i]=rho_NFW(R_grid[i]);
		psi_r[i]=-phi_NFW(R_grid[i]);
	}
#elif defined ISOCHRONE
	for (int i=0; i<NN; i++){
			rho_r[i]=rho_iso(R_grid[i]);
			psi_r[i]=-phi_iso(R_grid[i],1.);
	}
#endif

	/* gsl function: integrand (wrapper) */
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	gsl_function F;
	F.function = &Edd_DF::integrand_wrapper;
	F.params = this;

	/* initialize bar */
	printf("\n");
	ProgressBar bar(50,"=","notime");
	bar.init(NN);
	for (int i=1; i<NN; i++){
		bar.update(i+1);
		Ehere=psi_r[i];

		/* gsl integration */
		double err,res; int key=6;
		gsl_integration_qag (&F,0.,psi_r[i],1e-4,1e-3,5000,key,w,&res,&err);
		integral[i]=res;
	}
	integral[0]=integral[1]+(psi_r[0]-psi_r[1])/(psi_r[2]-psi_r[1])*(integral[2]-integral[1]);

	/* finalize bar */
	bar.fillSpace("\n..finishing DF(J) computation!!\n\n");

	gsl_integration_workspace_free (w);
}

/*
 *  integrand: returns the integrand of the formula, as a function of E.
 *  The derivative drho/dpsi is linearly approximated
 */
double Edd_DF::integrand(double psi){
	//return 1./sqrt(Ehere-psi)*dinterp(psi_r,rho_r,NN,psi);
	return 1./sqrt(Ehere-psi)*ddinterp(psi_r,rho_r,NN,psi);
}

/*
 *  dinterp: gsl approximation of first derivative
 *
 */
double Edd_DF::dinterp(double *xm,double *ym,int np,double x){ // linear interp for derivative
	int top,bot;
	if (topbottom(xm,np,x,&bot,&top)==0) return 0.;
	if (x<xm[0] || x>xm[np-1]) return 0;

	gsl_interp_init(interpolation, xm, ym, np);
	return gsl_interp_eval_deriv( interpolation, xm, ym, x, accelerator);
}

/*
 *  ddinterp: gsl approximation of second derivative
 *
 */
double Edd_DF::ddinterp(double *xm, double *ym, int np, double x){
	int top,bot;
	if (topbottom(xm,np,x,&bot,&top)==0) return 0.;
	if (bot==0) {bot++; top++;}
	if (top==np-1) {bot--; top--;}
	if (x<xm[0] || x>xm[np-1]) return 0;

	gsl_interp_init(interpolation, xm, ym, np);
	return gsl_interp_eval_deriv2 ( interpolation, xm, ym, x, accelerator );
}

/*
 *  integrate: calls the romberg (or trapezoidal) integration
 */
double Edd_DF::integrate(double E, int it){
	//double st;
	//return trapzd(0,E,&st,&it,1);

	return romberg(E,it);
}

/*
 *  trapzd: integrates integrand(double) from a to b in n intervals
 */
double Edd_DF::trapzd(double a, double b,double *s,int *it, int n){
	double x,sum,del;
	int j;
	if (n == 0) {
		*s=0.5*(b-a)*(integrand(a)+integrand(b));
		(*it)=1;
	} else {
		del=(b-a)/(*it);
		x=a+0.5*del; sum=0;
		for (j=0;j<(*it);j++,x+=del) sum += integrand(x);
		*s=0.5*(*s+(b-a)*sum/(*it));
		(*it)*=2;
	}
	return *s;
}

/*
 *  romberg: Romberg's rule for integration. Uses a combination of
 *  trapezoidal integrations.
 */
double Edd_DF::romberg(double E, int it){
	double st; int n=1;
	double int2,int1=trapzd(0,E,&st,&it,1);
	double res=1e6,tmp=0.,eps=0.00001;

	while(fabs(res-tmp)/res>eps){
		n++; tmp=int1;

		it=(int) it/2;
		int2=trapzd(0,E,&st,&it,1);
		res=( pow(4.,n-1)*int2 - int1 ) / (pow(4.,n-1) - 1.);

		if (n%10==0) eps*=10.;
	}

	return res;
}


/*
 *  edd_inversion: computes the DF as a function of E (input).
 *  Performs the last step of the formula, i.e., the derivative in E.
 *  (I use linear interpolation to estimate dg(E)/dE)
 */
double Edd_DF::edd_inversion(double E){
	return 1./(sqrt(8.)*M_PI*M_PI)*dinterp(psi_r,integral,NN,E);
}



/*
 * 	======================================
 *  GSL wrapper
 *  ======================================
 */
double Edd_DF::integrand_wrapper(double x, void *params)
    {
        return static_cast<Edd_DF*>(params)->integrand(x);
    }
