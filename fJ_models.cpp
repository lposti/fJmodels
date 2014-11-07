/*
 *Compute density using DF with first Chandra Phi and then Ylm Phi
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "press.h"
#include "potLEG.h"
#include "leg_pot3.h"
#include "isopot.h"
#include "tables3.h"
#include "uv_orb3.h"
#include "moments.h"
#include "check_acts.h"
#include "isodf4.h"
#include "print.h"

int rep=0;

void test(void){
	double R0=.001,z0;
	for(int i=0; i<10; i++){
		z0=.001;
		double r,M,FR,Fz,FR1,Fz1;
		//double dR=.001*R0, dz=.001*z0;
		//FR=-.5*(isopot(R0+dR,z0)-isopot(R0-dR,z0))/dR;
		//Fz=-.5*(isopot(R0,z0+dz)-isopot(R0,z0-dz))/dz;
		FR=-dPhiR(R0,z0); Fz=-dPhiz(R0,z0);
		isoforce(R0,z0,&FR1,&Fz1);
		r=sqrt(R0*R0+z0*z0); M=r*r*sqrt(FR1*FR1+Fz1*Fz1);
		printf("%8.4f %5.3f (%f %f) (%f %f) (%f %f)\n",R0,M,
		       phileg(R0,z0),phileg(R0,z0)/isopot(R0,z0)-1,
		       FR,FR/FR1-1,Fz,Fz/Fz1-1);
		R0*=1.7;
	}
}
void save_phil(double **phil_old,int nr,int np){
	for(int i=0; i<nr; i++)
		for(int j=0; j<np; j++) phil_old[i][j]=phil[i][j];
}
void save_Pr(double **Pr_old, double **Pr2_old, int nr,int np){
	for(int i=0; i<nr; i++)
		for(int j=0; j<np; j++) {
			Pr_old[i][j]=Pr[i][j];
			Pr2_old[i][j]=Pr2[i][j];
		}
}
double merge_phil(double **phil_old,double **Pr_old, double **Pr2_old,
				  int nr,int np,double f){
	double fa=1+f,maxd=0;
	printf("--merge_phil-> before: phil[%d][0]=%f phil[%d][1]=%f \n",nr-1,phil[nr-1][0],nr-1,phil[nr-1][1]);
	for(int i=0; i<nr; i++)
		for(int j=0; j<np; j++){
//			fprintf(efile,"%d %d %g %g\n",i,j,phil[i][j],phil_old[i][j]);
			maxd=MAX(maxd,fabs((phil[i][j]-phil_old[i][j])/phil[i][0]));
			phil[i][j]=fa*phil[i][j]-f*phil_old[i][j];
			Pr[i][j]=fa*Pr[i][j]-f*Pr_old[i][j];
			Pr2[i][j]=fa*Pr2[i][j]-f*Pr2_old[i][j];
		}
	return maxd;
}
/*
 *
 * Checks whether the grid (Lz,E,Er) can be successfully computed,
 * if not it restores the phil, Pr and Pr2 values to the latest computed
 * (i.e., implicitly assumes f=0 in merge_phil).
 *
 */
void check_pot(char *fname,double Rmax, double **phil_old2, double **Pr_old2,
				double **Pr2_old2,int nr, int npoly){
	if (tab_acts(fname,Rmax)!=0){
		printf("\n ****** Cannot build (Lz,E,Er) grid. Retrying "
				"with latest computed phil. ****** \n\n");
		for(int i=0; i<nr; i++)
				for(int j=0; j<npoly; j++){
					phil[i][j]=phil_old2[i][j];
					Pr[i][j]=Pr_old2[i][j];
					Pr2[i][j]=Pr2_old2[i][j];
				}
		if (tab_acts(fname,Rmax)!=0) exit(1);	// after restoring, re-run tab_acts
	}
}

int main(int nargs,char **args){
	/*
	if(nargs<3){
		printf("Must enter value of restart and nstep\n"); return 0;
	}
	int restart; sscanf(args[1],"%d",&restart);
	int nstep; sscanf(args[2],"%d",&nstep);
	*/
	int restart,nstep;
	time_t startTime = clock();
	double b=1,q=1.;
	isopot_init(b,q);//compute Phi of flattened isochrone for isopot()
	double Rmax=100*b;
	nr=120;	ar = new double[nr]; //allocate storage for potent()
	ngauss=6; npoly=3;
	double **phil_old, **phil_old2, **Pr_old, **Pr2_old, **Pr_old2, **Pr2_old2;
	phil_old=dmatrix(nr,npoly); phil_old2=dmatrix(nr,npoly);
	Pr_old=dmatrix(nr,npoly); Pr2_old=dmatrix(nr,npoly);
	Pr_old2=dmatrix(nr,npoly); Pr2_old2=dmatrix(nr,npoly);
	setgrid(.03*b,Rmax);
	rhl=dmatrix(nr,npoly); phil=dmatrix(nr,npoly);
	Pr=dmatrix(nr,npoly); Pr2=dmatrix(nr,npoly);

	double dr=1,dphi=1,dz=1;
	/*
	if (scanf("%d %d \n",&restart,&nstep)==0){
		printf("$$$ Problem reading first line of input file\n");
		exit(1);
	}
	if (scanf("%lf %lf %lf \n",&dr,&dphi,&dz)==0){
		printf("$$$ Problem reading second line of input file\n");
		exit(2);
	}
	*/
	char nameInput[40]; strcpy(nameInput,args[1]); printf("%s 1\n",nameInput);
	FILE *inpf=fopen(nameInput,"r");
	if (!fscanf(inpf,"%d %d",&restart,&nstep)) {
		printf("ERROR IN FIRST LINE INPUT!");
		exit(1);
	}
	printf("%d %d \n",restart,nstep);
	if (!fscanf(inpf,"%lf %lf %lf",&dr,&dphi,&dz)) {
		printf("ERROR IN SECOND LINE INPUT!");
		exit(1);
	}
	printf("%f %f %f \n",dr,dphi,dz);
#ifdef HJGJ
	double dr_g=1,dphi_g=1,dz_g=1;
	if (!fscanf(inpf,"%lf %lf %lf",&dr_g,&dphi_g,&dz_g)) {
		printf("ERROR IN THIRD LINE INPUT!");
		exit(1);
	}
	printf("%f %f %f\n",dr_g,dphi_g,dz_g);
#endif
	fclose(inpf);
	
#ifdef HJGJ
	setdf(dr,dphi,dz,dr_g,dphi_g,dz_g);
#else
	setdf(dr,dphi,dz);
#endif

	/* initial potential */
	setMJ0(1.,1.);
	phil_ini=dmatrix(nr,npoly); Pr_ini=dmatrix(nr,npoly); Pr2_ini=dmatrix(nr,npoly);
	char base[30],fname[30],stuff[30];
        strcpy(base,"models/hernq_newInPot_");
	int kontrl=1;


	if(restart!=0){
#ifdef HJGJ
		sprintf(stuff,"%4.2f_%4.2f_%4.2f_%4.2f_%d.out",dphi,dz,dphi_g,dz_g,restart-2);
#else
		sprintf(stuff,"%4.2f_%4.2f_%4.2f_%d.out",dr,dphi,dz,restart-2);
#endif
		strcpy(fname,base);strcat(fname,stuff);
		printf("Reading dump %s ..",fname);
		potent(fname,&rho,1,0);//compute Ylm coeffs using rho_l from fname1
		save_phil(phil_old,nr,npoly); save_Pr(Pr_old,Pr2_old,nr,npoly);
		sprintf(stuff,"_%d.out",restart-1);
		strcpy(fname,base);strcat(fname,stuff);
		printf("Reading dump %s\n",fname);
		potent(fname,&rho,1,0);//compute Ylm coeffs using rho_l from fname2
		printf("Restarting: Max shift: %7.4f\n",merge_phil(phil_old,Pr_old,Pr2_old,nr,npoly,0.));
		kontrl=1;
	}else{
#ifdef HJGJ
		sprintf(stuff,"%4.2f_%4.2f_%4.2f_%4.2f_-1.out",dphi,dz,dphi_g,dz_g);
#else
		sprintf(stuff,"%4.2f_%4.2f_%4.2f_-1.out",dr,dphi,dz);
#endif
		strcpy(fname,base);strcat(fname,stuff);
		printf("Computing Ylm Phi from analytic rho\n");
#ifdef HERNQUIST
		potent(fname,&rho_Hern,0,0);//compute Ylm coeffs for flattened Hernquist
#elif defined NFW
		potent(fname,&rho_NFW,0,0);//compute Ylm coeffs for flattened NFW
#elif defined ISOTHERMAL
		potent(fname,&rho_isoth,0,0);//compute Ylm coeffs for flattened Hernquist
#elif defined ISOCHRONE
		potent(fname,&isoden,0,0);//compute Ylm coeffs for flattened isochrone
#elif defined JAFFE
		potent(fname,&rho_Hern,0,0);//use Hernquist for the moment
#else
		printf("\n [ERROR]: must choose either Isochrone, Hernquist or NFW model in isodf4.h");
		exit(1);
#endif
		/* save initial potential */
		save_phil(phil_ini,nr,npoly); save_Pr(Pr_ini,Pr2_ini,nr,npoly);
	}
	tabstuff();

#ifdef PRINTDFH
	openDFH();
#endif
#ifdef PRINTXVJ
	openXVJ();
#endif

	for(rep=restart; rep<restart+nstep; rep++){
#ifdef HJGJ
		sprintf(stuff,"%4.2f_%4.2f_%4.2f_%4.2f_%d.acts",dphi,dz,dphi_g,dz_g,rep);
#else
		sprintf(stuff,"%4.2f_%4.2f_%4.2f_%d.acts",dr,dphi,dz,rep);
#endif
		strcpy(fname,base);strcat(fname,stuff);
		kontrl=1;
		if(kontrl==1){
			printf("tabulating JR(Lz,E,Iu) for n=\n");
			if (rep>0) check_pot(fname,Rmax,phil_old2,Pr_old2,Pr2_old2,nr,npoly);
			else tab_acts(fname,Rmax);
		} else {
			printf("reading JR(Lz,E,Er)..\n");
			get_acts(fname);
		}
        kontrl=1;

#ifdef HJGJ
		sprintf(stuff,"%4.2f_%4.2f_%4.2f_%4.2f_%d.out",dphi,dz,dphi_g,dz_g,rep);
#else
		sprintf(stuff,"%4.2f_%4.2f_%4.2f_%d.out",dr,dphi,dz,rep);
#endif
		strcpy(fname,base);strcat(fname,stuff);
		save_phil(phil_old,nr,npoly); save_Pr(Pr_old,Pr2_old,nr,npoly);
		printf("Computing Ylm Phi from DF..\n");
		potent5(fname,&rho,0,1);//compute Ylm coeffs using DF
#ifdef PRINTDFH
		closeDFH();
#endif
#ifdef PRINTXVJ
		closeXVJ();
#endif
		printf(">> %f %f %f",Pr_old2[0][0],Pr[0][0],Pr2[0][0]);
		save_phil(phil_old2,nr,npoly); save_Pr(Pr_old2,Pr2_old2,nr,npoly);

		/* CHANGED f PARAMETER IN merge_phil: NOW SET TO f=.25, before f=.5 */
		printf(" done. Max shift: %8.5f\n",merge_phil(phil_old,Pr_old,Pr2_old,nr,npoly,0.25));
		printf("-- time of computation: %f\n",(clock() - startTime) / (double) CLOCKS_PER_SEC);
	}
}
