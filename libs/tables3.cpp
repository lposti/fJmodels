//routines to set up and use tables of actions versus E and Er
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "potLEG.h"
#include "uv_orb3.h"
#include "grid3.h"
#include "tables3.h"
#include "press.h"
#include "mongo.h"
#include "leg_pot2.h"

/* Debug */
#include "Rzuv.h"

#define PIH 1.570796327
#define TINY 1.e-10
#define SMALL 1.e-5
//FILE *efile; 

#define SMALLEST 10

double get_closed(double,double,double);

int gEfn(int nL){
	return MAX(SMALLEST,NGRID*(float)(NR-nL)/NR);//more E values at low Lz
}
int gzfn(int nE){
	return MAX(SMALLEST,nE);//rises to NGRID with E}
}
double **tria(int j,int k){//more entries early on, fewer later
	double **x = new double*[j];
	for(int n=0; n<j; n++){
		int m=MAX(SMALLEST,k*(float)(j-n)/(j-1));
		x[n] = dmatrix(m);
	}
	return x;
}
double **tetra(int j,int k){//number of entries grows down the line
	double **x = new double*[j];// x[n] points to a line of values
	for(int n=0; n<j; n++){
		int m=MAX(SMALLEST,k*(float)n/(j-1));
		x[n] = dmatrix(m);
	}
	return x;
}	
double ***tetra(int i,int j,int k){
	double ***x = new double**[i];// x[n] points to a triangle of values
	for(int n=0; n<i; n++){
		int m=MAX(SMALLEST,j*(float)(i-n)/i);
		x[n]=tetra(m,k);
	}
	return x;
}
void compress_tria(FILE *ofile,double **x,int j,int k){
	for(int n=0; n<j; n++){
		int m=MAX(SMALLEST,k*(float)(j-n)/(j-1));
		compress(ofile,x[n],m);
	}
}
void get_tria(FILE *ifile,double **x,int j,int k){
	for(int n=0; n<j; n++){
		int m=MAX(SMALLEST,k*(float)(j-n)/(j-1));
		get(ifile,x[n],m);
		for(int i=0; i<m; i++) if((fabs(x[n][i])>1e50) || (isnan(x[n][i])==1)){
			printf("NAN in get_tria\n"); exit(0);
		}
	}
}
void compress_tetra(FILE *ofile,double **x,int j,int k){
	for(int n=0; n<j; n++){
		int m=MAX(SMALLEST,k*(float)n/(j-1));
		for(int i=0; i<m; i++) if((fabs(x[n][i])>1e50) || isnan(x[n][i])){
            printf("NAN in compress_tetra\n"); exit(0);
		}
		compress(ofile,x[n],m);
	}
}
void compress_tetra(FILE *ofile,double ***x,int i,int j, int k){
	for(int n=0; n<i; n++){
		int m=MAX(SMALLEST,j*(float)(i-n)/i);
		compress_tetra(ofile,x[n],m,k);
	}
}
void get_tetra(FILE *ifile,double **x,int j,int k){
	for(int n=0; n<j; n++){
		int m=MAX(SMALLEST,k*(float)n/(j-1));
		get(ifile,x[n],m);
		for(int i=0; i<m; i++) if((fabs(x[n][i])>1e50) || (isnan(x[n][i])==1)){
			printf("NAN in get_tetra\n"); exit(0);
		}
			
	}
}
void get_tetra(FILE *ifile,double ***x,int i,int j, int k){
	for(int n=0; n<i; n++){
		int m=MAX(SMALLEST,j*(float)(i-n)/i);
		get_tetra(ifile,x[n],m,k);
	}
}
void tabstuff(void){//call only once: just allocates storage
	Rgrid = dmatrix(NR); Lzgrid = dmatrix(NR); 
	Vmax = dmatrix(NR); Ecgrid = dmatrix(NR);
	Egrid = dmatrix(NGRID); Dgrid = dmatrix(NGRID);
	kappagrid = dmatrix(NR); nugrid = dmatrix(NR); Omegacgrid = dmatrix(NR);
	Jr_LzEIu=tetra(NR,NGRID,NGRID); Ergrid=tetra(NR,NGRID,NGRID);
	Jz_LzEIv=tetra(NR,NGRID,NGRID);
	sgrid=tria(NR,NGRID); Ermax = tria(NR,NGRID);
	//efile=fopen("tab_err2.dat","w");
}
static double RcEfn(double *Ep,double R){//what vanishes at Rc for given E
	return .5*R*dPhiR(R)+PhiR(R)-(*Ep);
}
static double RcE(double E,double Rc){
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
static double Rcfn(double *Lp,double R){//what vanishes at Rc for given Lz
	double x[2]={R,0},f[2]; dPhi(x,f);
	return (*Lp)/pow(R,3)-f[0];
}
static double GetRc(double Lzsq,double R){
	double Ri,Ro,dfRi,dfRo,dfRt=Rcfn(&Lzsq,R);
	if(dfRt>0){//inside
		while(dfRt>0){
			Ri=R; dfRi=dfRt; R*=1.2; dfRt=Rcfn(&Lzsq,R);
		} Ro=R; dfRo=dfRt;
	}else{//outside
		while(dfRt<0){
			Ro=R; dfRo=dfRt; R*=.8; dfRt=Rcfn(&Lzsq,R);
		} Ri=R; dfRi=dfRt;
	}
	return zbrent(&Lzsq,&Rcfn,Ri,Ro,dfRi,dfRo,TINY,20);
}
int choose_delta(void){
	//float Ri=1,zi=1; plots(1,&Ri,&zi,0,50,0,50,"R",1,"z",1,-.9,10);
	int gE=gEfn(1);
	for(int nE=0; nE<gE; nE++){//loop over E
		double E=Ecgrid[1]+.5*pow(Vmax[1]*(nE+1)/(double)gE,2);
		Dgrid[nE]=get_closed(RcE(E,Rgrid[1]),.25*Lzgrid[1],E);
		//printf("%g %g %g\n",E,Dgrid[nE],RcE(E,Rgrid[1]));
		Egrid[nE]=E;
	}
	//grend(1);
	return gE;
}
double Deltafn(double E){
	if(E<Egrid[0]) return sqrt(MAX(.000001,Dgrid[0]));
	else if(E<Egrid[nDgrid-1]) return sqrt(MAX(.000001,linterp(Egrid,Dgrid,nDgrid,E)));
	else return sqrt(MAX(.000001,Dgrid[nDgrid-1]));
}
int tab_acts(char *fname,double Rmax){//call on each change of Phi
	double Lzmax=Rmax*0.98*sqrt(-2*PhiR(Rmax)), Rmin=1e-3*Rmax;
	double Lzmin=Rmin*sqrt(Rmin*dPhiR(Rmin)),R=.01*Rmax;
    if (isnan(Lzmin)==1 || isnan(Lzmax)==1) {
    	printf("-- NAN Lzgrid: PhiR=%f dPhiR=%f\n",PhiR(Rmax),dPhiR(Rmin));
    	return 1;
    }
	for(int i=0; i<NR; i++){
		Lzgrid[i]=Lzmin+i*(Lzmax-Lzmin)/(double)(NR-1);
		double Lzsq=pow(Lzgrid[i],2); R=Rgrid[i]=GetRc(Lzsq,R);
		double vc2=R*dPhiR(R);
		Ecgrid[i]=.5*vc2+PhiR(R);//Ecirc at R
		Vmax[i]=.98*sqrt(-2*Ecgrid[i]);//Range of Egrid at R
		getfreqs(R,kappagrid+i,nugrid+i,Omegacgrid+i);
        if (Ecgrid[i]>0.) {
        	printf("-- Ecgrid[%d]=%f>0 @ R=%f PhiR=%f dPhiR=%f\n",i,Ecgrid[i],R,PhiR(R),dPhiR(R));
        	return 1;
        }
	}
	printf("Lz: (%f,%f,%f,..,%f)\n",Lzgrid[0],Lzgrid[1],Lzgrid[2],Lzgrid[NR-1]);
	nDgrid=choose_delta();
	//for(int i=0;i<nDgrid;i++) printf("%f %f\n",Egrid[i],Dgrid[i]);
	//plots(nDgrid,Egrid,Dgrid,Egrid[0],Egrid[nDgrid-1],-6,10,"E",1,"\\gD\\u2",6,-.9,10);
	//grend(1); if(nDgrid>0) exit(0);
	for(int nL=0; nL<NR; nL++){//loop over R & thus Lz
		if((nL % 10)==0) printf("%d ",nL);
		double Lz=Lzgrid[nL], R=Rgrid[nL];
		int gE=gEfn(nL),count;
		for(int nE=0; nE<gE; nE++){//loop over E
			double E=Ecgrid[nL]+.5*pow(Vmax[nL]*(nE+1)/(double)gE,2);
			double Delta=Deltafn(E);
			sgrid[nL][nE]=(double)(nE+1)/(double)gE;
			R*=3;
			double Phi0=PhiR(R), vc2=pow(Lz/R,2);
			double speed=sqrt(MAX(1e-6,2*(E-Phi0)-vc2));
            double utry=asinh(R/Delta);
			double x[2]={R,0}, p[2]={0,speed};
			uv_orb uvorb(Delta,Lz,Phi0,x,p);
			count=0;

            while(fabs(utry-uvorb.umid)>SMALL){//Find bottom of well U(u)
				x[0]=R=MAX(ar[0],Delta*sinh(uvorb.umid)); Phi0=PhiR(R);
				vc2=pow(Lz/R,2);
				p[1]=sqrt(MAX(1e-6,2*(E-Phi0)-vc2));
				utry=uvorb.umid;
                uvorb.reset(Delta,Lz,Phi0,x,p);
                /* Debug */

                count++;

                /*
                if (count>20) {
                	FILE *fU2=fopen("Upot2.dat","w");
                	fprintf(fU2,"# %f %f\n",uvorb.umid,uvorb.dU(uvorb.umid));
                	for (int nn=0;nn<1000;nn++){

                		double uhere=uvorb.umid-1.+nn*2./999.;
                		double xx[2]={Delta*sinh(uhere),0.};
                		double pp[2]={0.,sqrt(MAX(1e-6,2*(E-PhiR(xx[0]))-pow(Lz/xx[0],2)))};
                		pt X(Delta,xx,pp);
                		double ff[2],dphiu_here,dpuh=uvorb.dPhiu(uhere,ff);
                		dphiu_here=Delta*(ff[0]*cosh(uhere)*X.sv+ff[1]*sinh(uhere)*X.cv);
                		fprintf(fU2,"%f %f %f %f %f \n",uhere,uvorb.duturnfn(uhere),uvorb.uturnfn(uhere),dphiu_here,dpuh);
                	}
                	fclose(fU2); printf("nL=%d nE=%d \n",nL,nE);
                	utry=uvorb.umid;
                }
                */
			}
			speed=p[1];
			//Now fire off at each angle to plane
			int gz=gzfn(nE);
//#pragma omp parallel for

			for(int nz=0; nz<gz; nz++){
				double theta=nz*PIH/(double)(gz-1);
				double p[2]={speed*sin(theta),speed*cos(theta)};
				uv_orb uvorb(Delta,Lz,Phi0,x,p);
				Jr_LzEIu[nL][nE][nz]=uvorb.Ju();
				Jz_LzEIv[nL][nE][nz]=uvorb.Jv();
				if (uvorb.Er<-1.e-6) printf("Er=%f<0. !! \n",uvorb.Er);
				if (nz==gz-1 && uvorb.Er==0.) printf("Er=0. !!!! \n",uvorb.Er);
				Ergrid[nL][nE][nz]=sqrt(MAX(0,uvorb.Er));
			}
			Ermax[nL][nE]=Ergrid[nL][nE][gz-1];
            if (Ermax[nL][nE]==0.) {
            	printf("!!! Ermax[%d][%d]=0\n",nL,nE);
            	//for (int nz=0; nz<gz; nz++) printf("!!! Jr=%f Jz=%f\n",Jr_LzEIu[nL][nE][nz],Jz_LzEIv[nL][nE][nz]);
            	return 1;
            }
			for(int nz=0; nz<gz; nz++) {
				Ergrid[nL][nE][nz]/=Ermax[nL][nE];
				if (isnan(Ergrid[nL][nE][nz]==1)) {
					printf(">>NAN Ergrid[%d][%d][%d]\n",nL,nE,nz);
					return 1;
				}
			}
		}
	}
	printf("\n");
	FILE *acts_file; acts_file=fopen(fname,"w");
	if(acts_file==NULL){
		printf("I cannot open %s\n",fname); exit(0);
	}
	fprintf(acts_file,"%d\n",nDgrid);
	compress(acts_file,Egrid,nDgrid); compress(acts_file,Dgrid,nDgrid);
	compress(acts_file,Lzgrid,NR); compress(acts_file,Ecgrid,NR);
	compress(acts_file,Rgrid,NR); compress(acts_file,Vmax,NR);
	compress_tria(acts_file,Ermax,NR,NGRID);
	compress_tria(acts_file,sgrid,NR,NGRID);
    printf("compressing Er ..");
	compress_tetra(acts_file,Ergrid,NR,NGRID,NGRID);
	printf("compressing Jr ..");
	compress_tetra(acts_file,Jr_LzEIu,NR,NGRID,NGRID);
	printf("compressing Jz ..");
	compress_tetra(acts_file,Jz_LzEIv,NR,NGRID,NGRID);
	fclose(acts_file);
	return 0;
}
/*void check_tabs(void){
	for(int nL=0;nL<NR;nL++){
		for(int nE=1;nE<gEfn(nL);nE++){
			if(Ermax[nL][nE]<Ermax[nL][nE-1])
				fprintf(efile,"max: %d %d %f %f\n",nL,nE,Ermax[nL][nE],Ermax[nL][nE-1]);
			for(int nz=1;nz<gzfn(nE);nz++){
				if(Ergrid[nL][nE][nz]<Ergrid[nL][nE][nz-1])
					fprintf(efile,"Er %d %d %d %f %f\n",nL,nE,nz,Ergrid[nL][nE][nz],Ergrid[nL][nE][nz-1]);
			}
		}
	}
}*/
void get_acts(char *fname){//read in previously tabulated values
	FILE *acts_file; acts_file=fopen(fname,"r");
	if(acts_file==NULL){
		printf("I cannot open %s\n",fname); exit(0);
	}
	fscanf(acts_file,"%d",&nDgrid);
	get(acts_file,Egrid,nDgrid); get(acts_file,Dgrid,nDgrid);
	get(acts_file,Lzgrid,NR); get(acts_file,Ecgrid,NR);
	get(acts_file,Rgrid,NR); get(acts_file,Vmax,NR);
	get_tria(acts_file,Ermax,NR,NGRID);
	get_tria(acts_file,sgrid,NR,NGRID);
	get_tetra(acts_file,Ergrid,NR,NGRID,NGRID);
	get_tetra(acts_file,Jr_LzEIu,NR,NGRID,NGRID);
	get_tetra(acts_file,Jz_LzEIv,NR,NGRID,NGRID);
	fclose(acts_file);
	for(int i=0; i<NR; i++){
		getfreqs(Rgrid[i],kappagrid+i,nugrid+i,Omegacgrid+i);
	}
	//check_tabs();
}
int int_acts(double Lzin,double E,double Er,double *Jr,double *Jz){
	int b,t,bb,bt,tb,tt;
	double Lz=fabs(Lzin); if(Lz<Lzgrid[0] || Lz>Lzgrid[NR-1]){
		return 0;
	}
	//if(Lz<Lzgrid[0] || Lz>=0) return 0;//always abort look-up
	topbottom(Lzgrid,NR,Lz,&b,&t);//get position in Lz grid
	double dx=(Lz-Lzgrid[b])/(Lzgrid[t]-Lzgrid[b]), dxb=1-dx;
	double Ec=dxb*Ecgrid[b]+dx*Ecgrid[t];//estimated Ec(Lz)
	double sp=sqrt(MAX(0,2*(E-Ec)));//random speed
	double spm=dxb*Vmax[b]+dx*Vmax[t];//estimated max speed
	double s=sp/spm;
	if(s<0 || s>1){
		if(s<-0.02 || s>1.02) {
//			fprintf(efile,"s %d %f %f %f\n",b,s,Vmax[b],Vmax[t]);
			return 0;
		}
		s<0? s=0 : s=1;
	}
	double jubb,jubt,jutb,jutt,jvbb,jvbt,jvtb,jvtt,bdy,bdyb,tdy,tdyb,Ermbb,Ermtb;
	if(s<sgrid[b][0]){
		bt=0; bdy=s/sgrid[b][0]; bdyb=1-bdy; jubb=0; jvbb=0; Ermbb=0;
	}else{
		if(topbottom(sgrid[b],gEfn(b),s,&bb,&bt)==0) exit(0);
		bdy=(s-sgrid[b][bb])/(sgrid[b][bt]-sgrid[b][bb]), bdyb=1-bdy;
		Ermbb=Ermax[b][bb];
	}
	if(s<sgrid[t][0]){
		tt=0; tdy=s/sgrid[t][0]; tdyb=1-tdy; jutb=0; jvtb=0; Ermtb=0;
	}else{
		if(topbottom(sgrid[t],gEfn(t),s,&tb,&tt)==0) exit(0);
		tdy=(s-sgrid[t][tb])/(sgrid[t][tt]-sgrid[t][tb]); tdyb=1-tdy;
		Ermtb=Ermax[t][tb];
	}
	double Ermaxb=dxb*(bdyb*Ermbb+bdy*Ermax[b][bt])+dx*(tdyb*Ermtb+tdy*Ermax[t][tt]);
	double h=sqrt(MAX(0,Er))/MAX(1.e-20,Ermaxb);
	if(h<0 || h>1) {
		if(h<-.02 || h>1.02){
			if(s==1 && h>1){//we must be at or past Ve so kill DF
				*Jr=1e5; *Jz=0; return 1;
			}else{
				return 0;
			}
		}
		h<0? h=0 : h=1;
	}
	if(bt>0) jubb=linterp(Ergrid[b][bb],Jr_LzEIu[b][bb],gzfn(bb),h);
	jubt=linterp(Ergrid[b][bt],Jr_LzEIu[b][bt],gzfn(bt),h);
	if(tt>0) jutb=linterp(Ergrid[t][tb],Jr_LzEIu[t][tb],gzfn(tb),h);
	jutt=linterp(Ergrid[t][tt],Jr_LzEIu[t][tt],gzfn(tt),h);
	*Jr=dxb*(bdyb*jubb+bdy*jubt)+dx*(tdyb*jutb+tdy*jutt);//combine estimates of 4 cols
	if(bt>0) jvbb=linterp(Ergrid[b][bb],Jz_LzEIv[b][bb],gzfn(bb),h);
	jvbt=linterp(Ergrid[b][bt],Jz_LzEIv[b][bt],gzfn(bt),h);
	if(tt>0) jvtb=linterp(Ergrid[t][tb],Jz_LzEIv[t][tb],gzfn(tb),h);
	jvtt=linterp(Ergrid[t][tt],Jz_LzEIv[t][tt],gzfn(tt),h);
	*Jz=dxb*(bdyb*jvbb+bdy*jvbt)+dx*(tdyb*jvtb+tdy*jvtt);//combine estimates of 4 cols
	return 1;
}
