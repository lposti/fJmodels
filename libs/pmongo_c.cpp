#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN(A,B) ((A)<(B)?(A):(B))
#define MAX(A,B) ((A)>(B)?(A):(B))
#define ABS(A) ((A)<0?(-A):(A))
#define SIGN(A,B) ((B)>0?(A):(-A))
#define MOD(A,B) ((A)-(B)*((int)(A)/(B)))

static FILE *fp,*fpd;
static char first=1,isopen;
float *scratch;
double *dscratch;

void compress(FILE *tmpf,float *x,int npt){

	int i,j,k,lnx,rem,i1 = 90,i2 = 8100,i3 = 729000;
	float at2t23 = 8388608., big=1.e28,small=1.e-28,aln2;
	char m[6];
	m[5]='\0';

	aln2 = log(2.);
	for(k=0 ; k<npt ; k++){
		if (x[k]==0.){
			m[4] = 80;
			m[3] = 0;
			m[2] = 0;
			m[1] = 0;
			m[0] = 0;
		}else{
			if (ABS(x[k])>big) x[k]= SIGN(big, x[k]);
			if (ABS(x[k])<small) x[k] = SIGN(small, x[k]);
			if (ABS(x[k])>=1.)
				lnx = (int)(log(ABS(x[k]))/aln2) + 129;
			else
				lnx = (int)(log(ABS(x[k]))/aln2) + 128;
			j = lnx/i1;
			m[3] = lnx - j*i1;
			lnx = lnx - 129;
			if (lnx>0)
				rem = (int)((ABS(x[k])/pow(2.,(double)(lnx)) - 1.)*at2t23);
			else
				rem = (int)((ABS(x[k])*pow(2.,(double)(ABS(lnx))) - 1.)*at2t23);
			m[4] = rem/i3;
			rem = rem - m[4]*i3;
			m[2] = rem/i2;
			rem = rem - m[2]*i2;
			m[1] = rem/i1;
			m[0] = rem - m[1]*i1;
			m[4] = m[4] + j*12;
			if (x[k]<0.) m[4] = m[4] + 40;
		}
		for (i=0;i<5;i++){
			m[i] = m[i] + 33;
			if (m[i]>=94) m[i] = m[i] + 1;
			if (m[i]>=123) m[i] = m[i] + 1;
		}
		fprintf(tmpf,"%s",m);
		if(15*((k+1)/15)==k+1) fprintf(tmpf,"\n");
	}
	if(15*((npt)/15)!=npt) fprintf(tmpf,"\n");
}

void compres2(FILE *tmpf,float *array,int nx,int ny)
{
	int i,j;
	scratch=(float*)calloc(nx,sizeof(float));
	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++) scratch[j]=array[j+nx*i];
		compress(tmpf,scratch,nx);}
	free(scratch);
	scratch=NULL;
}
void compress(FILE *tmpf,double *x,int npt){

	int i,j,k,lnx,rem,i1 = 90,i2 = 8100,i3 = 729000;
	float at2t23 = 8388608., big=1.e28,small=1.e-28,aln2;
	char m[6];
	m[5]='\0';

	aln2 = log(2.);
	for(k=0 ; k<npt ; k++){
		if (x[k]==0.){
			m[4] = 80;
			m[3] = 0;
			m[2] = 0;
			m[1] = 0;
			m[0] = 0;
		}else{
			if (ABS(x[k])>big) x[k]= SIGN(big, x[k]);
			if (ABS(x[k])<small) x[k] = SIGN(small, x[k]);
			if (ABS(x[k])>=1.)
				lnx = (int)(log(ABS(x[k]))/aln2) + 129;
			else
				lnx = (int)(log(ABS(x[k]))/aln2) + 128;
			j = lnx/i1;
			m[3] = lnx - j*i1;
			lnx = lnx - 129;
			if (lnx>0)
				rem = (int)((ABS(x[k])/pow(2.,(double)(lnx)) - 1.)*at2t23);
			else
				rem = (int)((ABS(x[k])*pow(2.,(double)(ABS(lnx))) - 1.)*at2t23);
			m[4] = rem/i3;
			rem = rem - m[4]*i3;
			m[2] = rem/i2;
			rem = rem - m[2]*i2;
			m[1] = rem/i1;
			m[0] = rem - m[1]*i1;
			m[4] = m[4] + j*12;
			if (x[k]<0.) m[4] = m[4] + 40;
		}
		for (i=0;i<5;i++){
			m[i] = m[i] + 33;
			if (m[i]>=94) m[i] = m[i] + 1;
			if (m[i]>=123) m[i] = m[i] + 1;
		}
		fprintf(tmpf,"%s",m);
		if(15*((k+1)/15)==k+1) fprintf(tmpf,"\n");
	}
	if(15*((npt)/15)!=npt) fprintf(tmpf,"\n");
}

void get(FILE *tmpf,float *xg,int npt)
{
	float expon,at2t23= 8388608.;
	long int i1=90,i,j,k,line,nline;
	unsigned char lin[80],m[5];
	unsigned int np,kp;

	nline=npt/15;
	if(15*nline!=npt) nline=nline+1;
	for(line=0;line<nline;line++){
		fscanf(tmpf,"%s",lin);
		np=strlen((char*)lin)/5;
		if((5*np)!=strlen((char*)lin)) np=np+1;
		for(kp=0;kp<np;kp++){
			k=line*15+kp;
			for(i = 0;i<5;i++){
				m[i] = lin[5*kp+i];
				if (m[i]>=124) m[i] = m[i] - 1;
				if (m[i]>=95) m[i] = m[i] - 1;
				m[i] = m[i] - 33;}
			if (m[4]==80)
				xg[k] = 0.;
			else{
				if (m[4]>=40){
					m[4] = m[4] - 40;
					xg[k] = -1.;}
				else xg[k] = 1.;
				j = m[4]/12;
				m[4] = m[4] - j*12;
				expon=j*i1+m[3]-129;
				xg[k]=xg[k]*((((m[4]*i1+m[2])*i1+m[1])*i1
					      +m[0])/at2t23+ 1.);
				if (expon>0.)
					xg[k]=xg[k]*pow(2.,expon);
				else
					xg[k]=xg[k]/pow(2.,-expon);
			}
		}
	}
}

void get2(FILE *tmpf,float *array,int nx,int ny)
{
	int i,j;
	scratch=(float*)calloc(nx,sizeof(float));
	for(i=0;i<ny;i++){
		get(tmpf,scratch,nx);
		for(j=0;j<nx;j++) array[j+nx*i]=scratch[j];}
	free(scratch);
	scratch=NULL;
}
void compres2(FILE *tmpf,double *array,int nx,int ny)
{
	int i,j;
	dscratch=(double*)calloc(nx,sizeof(double));
	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++) dscratch[j]=array[j+nx*i];
		compress(tmpf,dscratch,nx);}
	free(dscratch);
	dscratch=NULL;
}
void compres2(FILE *tmpf,double **array,int nx,int ny){
	dscratch=(double*)calloc(nx,sizeof(double));
	for(int j=0;j<ny;j++){
		for(int i=0;i<nx;i++) dscratch[i]=array[i][j];
		compress(tmpf,dscratch,nx);
	}
	free(dscratch);
}
void compres2(FILE *tmpf,float **array,int nx,int ny){
	scratch=(float*)calloc(nx,sizeof(float));
	for(int j=0;j<ny;j++){
		for(int i=0;i<nx;i++) scratch[i]=array[i][j];
		compress(tmpf,scratch,nx);
	}
	free(scratch);
}
void compress(FILE *tmpf,double **array,int nx,int ny){//how it should be done
	for(int i=0; i<nx;i++) compress(tmpf,array[i],ny);
}
void compress(FILE *tmpf,double ***array,int nx,int ny,int nz){//might be compres3
	for(int i=0; i<nx; i++) compress(tmpf,array[i],ny,nz);
}
void get(FILE *tmpf,double *xg,int npt)
{
	float expon,at2t23= 8388608.;
	long int i1=90,i,j,k,line,nline;
	unsigned char lin[80],m[5];
	unsigned int np,kp;

	nline=npt/15;
	if(15*nline!=npt) nline=nline+1;
	for(line=0;line<nline;line++){
		fscanf(tmpf,"%s",lin);
		np=strlen((char*)lin)/5;
		if((5*np)!=strlen((char*)lin)) np=np+1;
		for(kp=0;kp<np;kp++){
			k=line*15+kp;
			for(i = 0;i<5;i++){
				m[i] = lin[5*kp+i];
				if (m[i]>=124) m[i] = m[i] - 1;
				if (m[i]>=95) m[i] = m[i] - 1;
				m[i] = m[i] - 33;}
			if (m[4]==80)
				xg[k] = 0.;
			else{
				if (m[4]>=40){
					m[4] = m[4] - 40;
					xg[k] = -1.;}
				else xg[k] = 1.;
				j = m[4]/12;
				m[4] = m[4] - j*12;
				expon=j*i1+m[3]-129;
				xg[k]=xg[k]*((((m[4]*i1+m[2])*i1+m[1])*i1
					      +m[0])/at2t23+ 1.);
				if (expon>0.)
					xg[k]=xg[k]*pow(2.,expon);
				else
					xg[k]=xg[k]/pow(2.,-expon);
			}
		}
	}
}

void get2(FILE *tmpf,double *array,int nx,int ny)
{
	int i,j;
	dscratch=(double*)calloc(nx,sizeof(double));
	for(i=0;i<ny;i++){
		get(tmpf,dscratch,nx);
		for(j=0;j<nx;j++) array[j+nx*i]=dscratch[j];}
	free(dscratch);
	dscratch=NULL;
}
void get2(FILE *tmpf,double **array,int nx,int ny){
	dscratch=(double*)calloc(nx,sizeof(double));
	for(int j=0;j<ny;j++){
		get(tmpf,dscratch,nx);
		for(int i=0; i<nx; i++) array[i][j]=dscratch[i];
	}
	free(dscratch);
}
void get(FILE *tmpf,double **array,int nx,int ny){//as it should be done
	for(int i=0; i<nx; i++) get(tmpf,array[i],ny);
}
void get(FILE *tmpf,double ***array,int nx,int ny,int nz){//might be get3
	for(int i=0; i<nx; i++) get(tmpf,array[i],ny,nz);
}
void choose(int kontrl,float *dx,float *dnx,float *dy,float *dny){
	int nx,ny,kch;

	if(kontrl==0){
		printf("Enter 1 1 for both scales log, 0 1 for log-y etc. ");
	    scanf("%i %i",&nx,&ny);
	    printf("Happy to let MONGO choose ticks (1)?");
	    scanf("%i",&kch);
	}else{
		fscanf(fpd,"%i %i",&nx,&ny);
        fscanf(fpd,"%i",&kch);
	}
	if(kch==1){
		*dx=-(float)nx;
		*dnx=-(float)nx;
		*dy=-(float)ny;
		*dny=-(float)ny;
	}else{
		if(kontrl==0){
			printf("Enter dx,dnx,dy,dny");
		    scanf("%g %g %g %g",dx,dnx,dy,dny);
		}else{
			fscanf(fpd,"%g %g %g %g",dx,dnx,dy,dny);
		}
		if(nx!=0)  *dx=-*dx;
		if(ny!=0)  *dy=-*dy;
	}
}

void grend(int nvec){
	fprintf(fp,"grend   \n");
	if(nvec==1) rewind(fpd);
}

void setltype(int m){
	fprintf(fp,"setltype\n");
	fprintf(fp,"%i\n",m);
}

void setlweight(int m){
	fprintf(fp,"setlweig\n");
	fprintf(fp,"%i\n",m);
}

void ctype(const char* color){
	fprintf(fp,"ctype   \n");
	fprintf(fp,"%s\n",color);
}

void setcolour(const char* color){
	fprintf(fp,"ctype   \n");
	fprintf(fp,"%s\n",color);
}

void relocate(float x,float y){
	fprintf(fp,"relocate\n");
	fprintf(fp,"%g %g\n",x,y);
}

void relocate(float *x){
	fprintf(fp,"relocate\n");
	fprintf(fp,"%g %g\n",x[0],x[1]);
}

void label(int m,const char *lab){
	fprintf(fp,"label   \n");
	fprintf(fp,"%i\n",m);
	fprintf(fp,"%s\n",lab);
}

void putlabel(int loc,int m,const char *lab){
/* Puts a label according to loc:
	7  8  9
	4  5  6
	1  2  3
 Left, center or right in x and y */
	fprintf(fp,"putlabel\n");
	fprintf(fp,"%i\n",loc);
	fprintf(fp,"%i\n",m);
	fprintf(fp,"%s\n",lab);
}

void setlim(float x1,float y1,float x2,float y2){
	fprintf(fp,"setlim  \n");
	fprintf(fp,"%g %g %g %g\n",x1,y1,x2,y2);
}

void window(int i,int j,int k){
	fprintf(fp,"window  \n");
	fprintf(fp,"%i %i %i\n",i,j,k);
}

void setangle(float a){
	fprintf(fp,"setangle\n");
	fprintf(fp,"%g\n",a);
}

void point(float a){//4 for arrows
	fprintf(fp,"point   \n");
	fprintf(fp,"%8.2f\n",a);
}

void points(float a,int ns,float *x,float *y,int n){
	int i,np,nr,nstart;
	nr=n/512;
	for(i=0;i<nr;i++){
		fprintf(fp,"points  \n");
		fprintf(fp,"%i\n",512);
		fprintf(fp,"%g %i\n",a,ns);
		nstart=512*i;
		compress(fp,&x[nstart],512);
		compress(fp,&y[nstart],512);
	}
	np=n-nr*512;
	if(np>0){
		fprintf(fp,"points  \n");
		fprintf(fp,"%i\n",np);
		fprintf(fp,"%g %i\n",a,ns);
		nstart=512*nr;
		compress(fp,&x[nstart],np);
		compress(fp,&y[nstart],np);
	}
}
void points(double a,int ns,double *x,double *y,int n){
	int i,np,nr,nstart;
	nr=n/512;
	for(i=0;i<nr;i++){
		fprintf(fp,"points  \n");
		fprintf(fp,"%i\n",512);
		fprintf(fp,"%g %i\n",a,ns);
		nstart=512*i;
		compress(fp,&x[nstart],512);
		compress(fp,&y[nstart],512);
	}
	np=n-nr*512;
	if(np>0){
		fprintf(fp,"points  \n");
		fprintf(fp,"%i\n",np);
		fprintf(fp,"%g %i\n",a,ns);
		nstart=512*nr;
		compress(fp,&x[nstart],np);
		compress(fp,&y[nstart],np);
	}
}

void setexpand(float a){
	fprintf(fp,"setexpan\n");
	fprintf(fp,"%g\n",a);
}

void draw(float a,float b){
	fprintf(fp,"draw    \n");
	fprintf(fp,"%g %g\n",a,b);
}

void draw(float *a){
	fprintf(fp,"draw    \n");
	fprintf(fp,"%g %g\n",a[0],a[1]);
}

void connect(float *x,float *y,int n){
	fprintf(fp,"connect \n");
	fprintf(fp,"%i\n",n);
	compress(fp,x,n);
	compress(fp,y,n);
}

void connect(double *x,double *y,int n){
	fprintf(fp,"connect \n");
	fprintf(fp,"%i\n",n);
	compress(fp,x,n);
	compress(fp,y,n);
}

void mgobox(int m,int n,int n1,int n2){
	fprintf(fp,"mgobox  \n");
	fprintf(fp,"%i %i %i %i\n",m,n,n1,n2);
}

void xlabel(int m,const char *lab){
	fprintf(fp,"xlabel  \n");
	fprintf(fp,"%i\n",m);
	fprintf(fp,"%s\n",lab);
}

void ylabel(int m,const char *lab){
	fprintf(fp,"ylabel  \n");
	fprintf(fp,"%i\n",m);
	fprintf(fp,"%s\n",lab);
}

void errorbar(int location,float *x,float *y,float *err,int nxy){
/*
c location =1,2,3,4, for +x,+y,-x,-y
*/
	fprintf(fp,"errorbar\n");
	fprintf(fp,"%i %i\n",location,nxy);
	compress(fp,x,nxy);
	compress(fp,y,nxy);
	compress(fp,err,nxy);
}

void errorbar(int location,double *x,double *y,double *err,int nxy){
/*
c location =1,2,3,4, for +x,+y,-x,-y
*/
	fprintf(fp,"errorbar\n");
	fprintf(fp,"%i %i\n",location,nxy);
	compress(fp,x,nxy);
	compress(fp,y,nxy);
	compress(fp,err,nxy);
}

void contour(float *array,int nx,int ny,float *h,int ncont){
	fprintf(fp,"contour \n");
	fprintf(fp,"%i %i %i\n",nx,ny,ncont);
	compress(fp,h,ncont);
	compres2(fp,array,nx,ny);
}
void contour(double **array,int nx,int ny,double *h,int ncont){
	fprintf(fp,"contour \n");
	fprintf(fp,"%i %i %i\n",nx,ny,ncont);
	compress(fp,h,ncont);
	compres2(fp,array,nx,ny);
}
void contour(float **array,int nx,int ny,float *h,int ncont){
	fprintf(fp,"contour \n");
	fprintf(fp,"%i %i %i\n",nx,ny,ncont);
	compress(fp,h,ncont);
	compres2(fp,array,nx,ny);
}
void scontour(float *z,int nx,int ny,float *level,int nlevel,float theta,float phi, float chi,int ifaxis){
	fprintf(fp,"scontour\n");
	fprintf(fp,"%i %i %i\n",nx,ny,nlevel);
	compress(fp,level,nlevel);
	compres2(fp,z,nx,ny);
	fprintf(fp,"%f %f %f %i\n",theta,phi,chi,ifaxis);
}

void histogram(float *x,float *num,int npt){
	fprintf(fp,"histogra\n");
	fprintf(fp,"%i\n",npt);
	compress(fp,x,npt);
	compress(fp,num,npt);
}

void histogram(double *x,double *num,int npt){
	fprintf(fp,"histogra\n");
	fprintf(fp,"%i\n",npt);
	compress(fp,x,npt);
	compress(fp,num,npt);
}

void plt3d(float *a,int m,int n,float alt,float az,float zfac,float zoff){
	fprintf(fp,"plt3d   \n");
	fprintf(fp,"%i %i\n",m,n);
	compres2(fp,a,m,n);
	fprintf(fp,"%g %g %g %g\n",alt,az,zfac,zoff);
}

void plt3d(double **a,int m,int n,double alt,double az,double zfac,double zoff){
	fprintf(fp,"plt3d   \n");
	fprintf(fp,"%i %i\n",m,n);
	compres2(fp,a,m,n);
	fprintf(fp,"%g %g %g %g\n",alt,az,zfac,zoff);
}

void label3d(const char *xlab,int nx,const char *ylab,int ny,const char *zlab,int nz){
	fprintf(fp,"label3d \n");
	fprintf(fp,"%i %i %i\n",nx,ny,nz);
	fprintf(fp,"%s\n",xlab);
	fprintf(fp,"%s\n",ylab);
	fprintf(fp,"%s\n",zlab);
}

/*void arrows(float *x,int np,float scale){//x[0]=x,x[1]=y,x[2]=|v|,x[3]=theta
	int i=0; float rad=0.0174533;
	while(i<np){
		relocate(x[i]-x[i+2]*scale*cos(rad*x[i+3]),x[i+1]-x[i+2]*scale*sin(rad*x[i+3]));
		draw(x[i],x[i+1]);
		setangle(x[i+3]);
		point(54.+.1*x[i+2]);
		i+=4;
	}
}*/
void arrows(float *x,int npoint,float scale){//x[0]=x,x[1]=y,x[2]=|v|,x[3]=theta
	fprintf(fp,"arrows  \n");
	fprintf(fp,"%d %f\n",npoint,scale);
	compress(fp,x,4*npoint);
}

void greyscal(float *array,int nx,int ny,float white,float black){
	fprintf(fp,"greyscal\n");
	fprintf(fp,"%i %i %g %g\n",nx,ny,white,black);
	compres2(fp,array,nx,ny);
}
void greyscal(double **array,int nx,int ny,double white,double black){
	fprintf(fp,"greyscal\n");
	fprintf(fp,"%i %i %lg %lg\n",nx,ny,white,black);
	compres2(fp,array,nx,ny);
}

void coloursc(float* array,int nx,int ny,float cold,float hot){
	fprintf(fp,"coloursc\n");
	fprintf(fp,"%i %i %g %g\n",nx,ny,cold,hot);
	compres2(fp,array,nx,ny);
}
void coloursc(double** array,int nx,int ny,double cold,double hot){
	fprintf(fp,"coloursc\n");
	fprintf(fp,"%i %i %lg %lg\n",nx,ny,cold,hot);
	compres2(fp,array,nx,ny);
}
void colourwh(float* array,int nx,int ny,float cold,float hot,float white){
	fprintf(fp,"colourwh\n");
	fprintf(fp,"%i %i %g %g %g\n",nx,ny,cold,hot,white);
	compres2(fp,array,nx,ny);
}
void colourwh(float** array,int nx,int ny,float cold,float hot,float white){
	fprintf(fp,"colourwh\n");
	fprintf(fp,"%i %i %g %g %g\n",nx,ny,cold,hot,white);
	compres2(fp,array,nx,ny);
}
void colourwh(double** array,int nx,int ny,double cold,double hot,double white){
	fprintf(fp,"colourwh\n");
	fprintf(fp,"%i %i %g %g %g\n",nx,ny,cold,hot,white);
	compres2(fp,array,nx,ny);
}
void colourwR(float* array,int nx,int ny,float cold,float hot,float white){
	fprintf(fp,"colourwR\n");
	fprintf(fp,"%i %i %g %g %g\n",nx,ny,cold,hot,white);
	compres2(fp,array,nx,ny);
}
void colourwR(float** array,int nx,int ny,double cold,double hot,double white){
    fprintf(fp,"colourwR\n");
    fprintf(fp,"%i %i %g %g %g\n",nx,ny,cold,hot,white);
    compres2(fp,array,nx,ny);
}
void colourwR(double** array,int nx,int ny,double cold,double hot,double white){
	fprintf(fp,"colourwR\n");
	fprintf(fp,"%i %i %g %g %g\n",nx,ny,cold,hot,white);
	compres2(fp,array,nx,ny);
}
void fillin(){
	fprintf(fp,"fillin  \n");
}

void rgbcolor(int i,int j,int k){
	fprintf(fp,"rgbcolor\n");
	fprintf(fp,"%i %i %i\n",i,j,k);
}

void shade(int delta,float *x,float *y,int np){
	fprintf(fp,"shade   \n");
	fprintf(fp,"%i %i",np,delta);
	compress(fp,x,np);
	compress(fp,y,np);
}
float plhist(float bottom,float *upper,float *n,int np){//plots dN/ds
	//bottom lowest x, n[i] # objects between bottom or upper[i-1]
	//& upper[i] returns total # of objects
	//  in calling pgm: j=0; while(diff>upper[j] && j<np-1) j++; n[j]+=1;
	float *errx = new float[np]; float *erry = new float[np];
	float *plotx = new float[np]; float *ploty = new float[np];
	plotx[0]=.5*(bottom+upper[0]);
	ploty[0]=n[0]/(upper[0]-bottom);
	errx[0]=upper[0]-plotx[0];
	erry[0]=sqrt(n[0])/(upper[0]-bottom);
	float total=n[0];
	for(int i=1;i<np;i++){
		total+=n[i];
		plotx[i]=.5*(upper[i-1]+upper[i]);
		ploty[i]=n[i]/(upper[i]-upper[i-1]);
		errx[i]=upper[i]-plotx[i];
		erry[i]=sqrt(n[i])/(upper[i]-upper[i-1]);
	}
	points(40.13,1,plotx,ploty,np);
	errorbar(1,plotx,ploty,errx,np); errorbar(3,plotx,ploty,errx,np);
	errorbar(2,plotx,ploty,erry,np); errorbar(4,plotx,ploty,erry,np);
	delete [] errx; delete [] erry; delete [] plotx; delete [] ploty;
	return total;
}
#define TINY 1e-20
float plloghist(float bottom,float *upper,float *n,int np){//plots log(dN/ds)
	//bottom lowest x, n[i] # objects between bottom or upper[i-1] & upper[i]
	float *errx = new float[np]; float *erry = new float[np];
	float *plotx = new float[np]; float *ploty = new float[np];
	plotx[0]=.5*(bottom+upper[0]);
	ploty[0]=log10(MAX(TINY,n[0])/(upper[0]-bottom));
	errx[0]=upper[0]-plotx[0];
	erry[0]=log10(1+1/sqrt(MAX(TINY,n[0])));
	float total=n[0];
	for(int i=1;i<np;i++){
		total+=n[i];
		plotx[i]=.5*(upper[i-1]+upper[i]);
		ploty[i]=log10(MAX(TINY,n[i])/(upper[i]-upper[i-1]));
		errx[i]=upper[i]-plotx[i];
		erry[i]=log10(1+1/sqrt(MAX(TINY,n[i])));
	}
	points(40.13,1,plotx,ploty,np);
	errorbar(1,plotx,ploty,errx,np); errorbar(2,plotx,ploty,erry,np);
	for(int i=0;i<np;i++)
		erry[i]=-log10(MAX(1.e-20,1-1/sqrt(MAX(TINY,n[i]))));
	errorbar(3,plotx,ploty,errx,np); errorbar(4,plotx,ploty,erry,np);
	delete [] errx; delete [] erry; delete [] plotx; delete [] ploty;
	return total;
}
#undef TINY
void plots(int nmax,float *x,float *y,float xmin0,float xmax0,float
	   ymin0,float ymax0,const char *xlab,int nxl,const char *ylab,int nyl,
	   float ex,int ikont){

	int i,j,k,kontrl,ksq,l,ndev,ngrend=0,npane,nwind;
	float dx,dnx,dy,dny,yox;
	char fname[30];

	kontrl=ikont;
	if(first){
		first=0;
		isopen=0;
		strcpy(fname,getenv("mdir"));
		strcat(fname,"\\mongplot.tmp");
		if((fp=fopen(fname,"w"))==NULL){
			printf("I cannot open %s\n",fname);
			exit(0);
		}
	}
	fprintf(fp,"plots   \n");
	fprintf(fp,"%i\n",nmax);
	compress(fp,x,nmax);
	compress(fp,y,nmax);
	k=MOD(kontrl,10);
	fprintf(fp,"%g %g %g %g %i %i %g %i\n",xmin0,xmax0,ymin0,ymax0,nxl,nyl,ex,k);
	if(nxl>0)  fprintf(fp,"%s\n",xlab);
	if(nyl>0)  fprintf(fp,"%s\n",ylab);
	if(ngrend==0)  nwind=0;
	if(kontrl==10){
		if(isopen==0){
			if((fpd=fopen("default.plt","r"))==NULL){
				printf("failed to open def\n");fflush(stdout);
				isopen=0;
				kontrl=0;
				printf("Number of WINDOW panes? ");
				scanf("%i",&nwind);
				fprintf(fp,"%i\n",nwind);
				if(nwind!=1){
					printf("Which pane should I be plotting in?");
					scanf("%i",&npane);
					fprintf(fp,"%i\n",npane);
				}
				printf("Should x- and y-axes have same scale (1)? (2 for specified y/x, 4 for square)");
				scanf("%i",&ksq);
				fprintf(fp,"%i\n",ksq);
				if(ksq==2){
					printf("Required value of y/x?");
					scanf("%g",&yox);
					fprintf(fp,"%g\n",yox);
				}
				if(nxl+nyl>0){
					choose(kontrl,&dx,&dnx,&dy,&dny);
					fprintf(fp,"%g %g %g %g\n",dx,dnx,dy,dny);
				}else{
					if(ndev!=1) printf("No box plotted because no labels given");
				}
			}else
				isopen=1;
		}
		if(isopen==1){
			fscanf(fpd,"%i",&nwind);
			fprintf(fp,"%i\n",nwind);
			if(nwind!=1){
				if(ABS(nwind)>2){
					fscanf(fpd,"%i %i %i %i",&i,&j,&k,&l);
					fprintf(fp,"%i %i %i %i\n",i,j,k,l);
				}else{
					fscanf(fpd,"%i",&npane);
					fprintf(fp,"%i\n",npane);
				}
			}
			fscanf(fpd,"%i",&ksq);
			fprintf(fp,"%i\n",ksq);
			if(ksq==2){
				fscanf(fpd,"%g",&yox);
				fprintf(fp,"%g\n",yox);
			}
			if(nxl+nyl>0){
				choose(kontrl,&dx,&dnx,&dy,&dny);
				fprintf(fp,"%g %g %g %g\n",dx,dnx,dy,dny);
			}else{
				if(ndev!=1)  printf("No box plotted because no labels given\n");
			}
		}
//		printf("leaving plots\n");fflush(stdout);
	}
}
void plots(int nmax,double *x,double *y,float xmin0,float xmax0,float
	   ymin0,float ymax0,const char *xlab,int nxl,const char *ylab,int nyl,
	   float ex,int ikont){

	int i,j,k,kontrl,ksq,l,ndev,ngrend=0,npane,nwind;
	float dx,dnx,dy,dny,yox;
	char fname[30];

	kontrl=ikont;
	if(first){
		first=0;
		isopen=0;
        /*
        * strcpy(fname,getenv("mdir"));
        * strcat(fname,"\\mongplot.tmp");
        */
        strcpy(fname,"mongplot.tmp");
		if((fp=fopen(fname,"w"))==NULL){
			printf("I cannot open %s\n",fname);
			exit(0);
		}
	}
	fprintf(fp,"plots   \n");
	fprintf(fp,"%i\n",nmax);
	compress(fp,x,nmax);
	compress(fp,y,nmax);
	k=MOD(kontrl,10);
	fprintf(fp,"%g %g %g %g %i %i %g %i\n",xmin0,xmax0,ymin0,ymax0,nxl,nyl,ex,k);
	if(nxl>0)  fprintf(fp,"%s\n",xlab);
	if(nyl>0)  fprintf(fp,"%s\n",ylab);
	if(ngrend==0)  nwind=0;
	if(kontrl==10){
		if(isopen==0){
			if((fpd=fopen("default.plt","r"))==NULL){
				printf("failed to open def\n");fflush(stdout);
				isopen=0;
				kontrl=0;
				printf("Number of WINDOW panes? ");
				scanf("%i",&nwind);
				fprintf(fp,"%i\n",nwind);
				if(nwind!=1){
					printf("Which pane should I be plotting in?");
					scanf("%i",&npane);
					fprintf(fp,"%i\n",npane);
				}
				printf("Should x- and y-axes have same scale (1)? (2 for specified y/x, 4 for square)");
				scanf("%i",&ksq);
				fprintf(fp,"%i\n",ksq);
				if(ksq==2){
					printf("Required value of y/x?");
					scanf("%g",&yox);
					fprintf(fp,"%g\n",yox);
				}
				if(nxl+nyl>0){
					choose(kontrl,&dx,&dnx,&dy,&dny);
					fprintf(fp,"%g %g %g %g\n",dx,dnx,dy,dny);
				}else{
					if(ndev!=1) printf("No box plotted because no labels given");
				}
			}else
				isopen=1;
		}
		if(isopen==1){
			fscanf(fpd,"%i",&nwind);
			fprintf(fp,"%i\n",nwind);
			if(nwind!=1){
				if(ABS(nwind)>4){
					fscanf(fpd,"%i %i %i %i",&i,&j,&k,&l);
					fprintf(fp,"%i %i %i %i\n",i,j,k,l);
				}else{
					fscanf(fpd,"%i",&npane);
					fprintf(fp,"%i\n",npane);
				}
			}
			fscanf(fpd,"%i",&ksq);
			fprintf(fp,"%i\n",ksq);
			if(ksq==2){
				fscanf(fpd,"%g",&yox);
				fprintf(fp,"%g\n",yox);
			}
			if(nxl+nyl>0){
				choose(kontrl,&dx,&dnx,&dy,&dny);
				fprintf(fp,"%g %g %g %g\n",dx,dnx,dy,dny);
			}else{
				if(ndev!=1)  printf("No box plotted because no labels given\n");
			}
		}
//		printf("leaving plots\n");fflush(stdout);
	}
}
