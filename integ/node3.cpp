#include <math.h>
#include "node3.h"

#define SMALL 1.e-7

#define TMPARRAY tmp[0]=c_value[0]; tmp[1]=c_value[1]; tmp[2]=c_value[2]; tmp[3]=c_value[3];\
				 tmp[4]=c_value[4]; tmp[5]=c_value[5]; tmp[6]=c_value[6]; tmp[7]=c_value[7]
#define TMP2ARRAY(XX) tmp2[0]=fix3[0].XX*diff.XX; tmp2[1]=fix3[1].XX*diff.XX; tmp2[2]=fix3[2].XX*diff.XX; tmp2[3]=fix3[3].XX*diff.XX;\
					  tmp2[4]=fix3[4].XX*diff.XX; tmp2[5]=fix3[5].XX*diff.XX; tmp2[6]=fix3[6].XX*diff.XX; tmp2[7]=fix3[7].XX*diff.XX
#define TMP3ARRAY(XX) tmp3[0]=fix3[0].XX*diff.XX; tmp3[1]=fix3[1].XX*diff.XX; tmp3[2]=fix3[2].XX*diff.XX; tmp3[3]=fix3[3].XX*diff.XX;\
					  tmp3[4]=fix3[4].XX*diff.XX; tmp3[5]=fix3[5].XX*diff.XX; tmp3[6]=fix3[6].XX*diff.XX; tmp3[7]=fix3[7].XX*diff.XX


intloc3 fix3[8]={{0,0,0},{0,1,0},{1,1,0},{1,0,0},{0,0,1},{0,1,1},{1,1,1},{1,0,1}};

double node3::GetSum(void){
	double sum=0.,tmp[8];
	TMPARRAY;
	for(int i=0;i<8;i++) sum+=tmp[i];
	return sum;
}
double node3::GetVariance(double * __restrict__ sum1){
	//double sum2=pow(c_value[0],2); *sum1=c_value[0];
	double sum=0.,sum2=0.,tmp[8];
	TMPARRAY;
	for(unsigned i=0;i<8;i++){
		sum+=tmp[i]; sum2+=tmp[i]*tmp[i];
	}
	*sum1=sum;
	return .125*(sum2-.125*pow(*sum1,2));
}
double node3::GetI(void){// evaluation of rho
	double sum=0.,tmp[8];
	TMPARRAY;
	for(int i=0;i<8;i++) sum+=tmp[i];
	return .5*(c_value[8]+.125*sum);
}
double node3::GetIyb(void){//for evaluation of rho*ybar 
	double sum=0,tmp[8],tmp2[8];
	TMPARRAY;
	TMP2ARRAY(y);
	for(int i=0;i<8;i++) sum+=tmp[i]*(loc.y+tmp2[i]);
	return .5*(c_value[8]*(loc.y+0.5*diff.y)+.125*sum);
}
double node3::GetIfy(double *x,double (*fL)(double)){//for evaluation of rho*ybar 
	double sum=0;
	for(int i=0;i<8;i++){
		double v=loc.y+fix3[i].y*diff.y;
		sum+=c_value[i]*v*(*fL)(x[0]*v);
	}
	double v=loc.y+0.5*diff.y;
	return .5*(c_value[8]*v*(*fL)(x[0]*v)+.125*sum);	// to change rotation sign, change sign of this expression
}														// see Binney (2014) eq 16
double node3::GetIxsq(void){//for evaluation of rho*xsqbar
	double sum=0,tmp[8],tmp2[8];
	TMPARRAY;
	TMP2ARRAY(x);
	for(int i=0;i<8;i++) sum+=tmp[i]*(loc.x+tmp2[i])*(loc.x+tmp2[i]);
	return .5*(c_value[8]*pow(loc.x+0.5*diff.x,2)+.125*sum);
}
double node3::GetIysq(void){//for evaluation of rho*ysqbar
	double sum=0,tmp[8],tmp2[8];
	TMPARRAY;
	TMP2ARRAY(y);
	for(int i=0;i<8;i++) sum+=tmp[i]*(loc.y+tmp2[i])*(loc.y+tmp2[i]);
	return .5*(c_value[8]*pow(loc.y+0.5*diff.y,2)+.125*sum);
}
double node3::GetIzsq(void){//for evaluation of rho*zsqbar
	double sum=0,tmp[8],tmp2[8];
	TMPARRAY;
	TMP2ARRAY(z);
	for(int i=0;i<8;i++) sum+=tmp[i]*(loc.z+tmp2[i])*(loc.z+tmp2[i]);
	return .5*(c_value[8]*pow(loc.z+0.5*diff.z,2)+.125*sum);
}
double node3::GetIxz(void){//for evaluation of rho*zsqbar
	double sum=0,tmp[8],tmp2[8],tmp3[8];
	TMPARRAY;
	TMP2ARRAY(x);
	TMP3ARRAY(z);
	for(int i=0;i<8;i++) sum+=tmp[i]*(loc.x+tmp2[i])*(loc.z+tmp3[i]);
	return .5*(c_value[8]*(loc.x+0.5*diff.x)*(loc.z+0.5*diff.z)+.125*sum);
}
void node3::GetCvals(double *vals){
	for(int i=0;i<9;i++) vals[i]=c_value[i];
}
void node3::GetVerts(location3 *verts){
	for(int i=0;i<8;i++){
		verts[i].x=loc.x+fix3[i].x*diff.x;
		verts[i].y=loc.y+fix3[i].y*diff.y;
		verts[i].z=loc.z+fix3[i].z*diff.z;
	}
	verts[8].x=loc.x+.5*diff.x; verts[8].y=loc.y+.5*diff.y; verts[8].z=loc.z+.5*diff.z;
}
void node3::set_values(double (*fn)(double,double,double)){
	location3 points[9];
	GetVerts(points);
	for(int i=0;i<9;i++) c_value[i]=(*fn)(points[i].x,points[i].y,points[i].z);
}
void node3::set_values(double *x,double (*fn)(double*,double,double,double)){
	location3 points[9];
	GetVerts(points);
	for(int i=0;i<9;i++) c_value[i]=(*fn)(x,points[i].x,points[i].y,points[i].z);
}
node3::node3(location3 L,int lev,location3 DIFF,double value0,double value1,double value2,double value3,
	     double value4,double value5,double value6,double value7,double value8){
	loc=L; level=lev;
	double k=pow(2,level);
	diff.x=DIFF.x/k; diff.y=DIFF.y/k; diff.z=DIFF.z/k;
	c_value[0]=value0; c_value[1]=value1; c_value[2]=value2; c_value[3]=value3;
	c_value[4]=value4; c_value[5]=value5; c_value[6]=value6; c_value[7]=value7;
	c_value[8]=value8;
}
node3::node3(double lx,double ly,double lz,int lev,location3 DIFF){
	loc.x=lx; loc.y=ly; loc.z=lz; level=lev;
	double k=pow(2,level);
	diff.x=DIFF.x/k; diff.y=DIFF.y/k; diff.z=DIFF.z/k;
}
node3::~node3(void){}
void eval_node3(double (*fn)(double,double,double),node3 parent,location3 DIFF,
		double *I,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for(int i=0;i<8;i++) m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from  
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax) *I+=I1*dx*dy*dz;
	else{
		eval_node3(fn,kiddie0,DIFF,I,lmax); eval_node3(fn,kiddie1,DIFF,I,lmax);
		eval_node3(fn,kiddie2,DIFF,I,lmax); eval_node3(fn,kiddie3,DIFF,I,lmax);
		eval_node3(fn,kiddie4,DIFF,I,lmax); eval_node3(fn,kiddie5,DIFF,I,lmax);
		eval_node3(fn,kiddie6,DIFF,I,lmax); eval_node3(fn,kiddie7,DIFF,I,lmax);
	}
}
void eval_node3(double *Rz,double (*fn)(double*,double,double,double),node3 parent,location3 DIFF,
		double *I,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for(int i=0;i<8;i++) m[i]=(*fn)(Rz,BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(Rz,BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from  
	f1=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(Rz,BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(Rz,BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(Rz,BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(Rz,BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax) *I+=I1*dx*dy*dz;
	else{
		eval_node3(Rz,fn,kiddie0,DIFF,I,lmax); eval_node3(Rz,fn,kiddie1,DIFF,I,lmax);
		eval_node3(Rz,fn,kiddie2,DIFF,I,lmax); eval_node3(Rz,fn,kiddie3,DIFF,I,lmax);
		eval_node3(Rz,fn,kiddie4,DIFF,I,lmax); eval_node3(Rz,fn,kiddie5,DIFF,I,lmax);
		eval_node3(Rz,fn,kiddie6,DIFF,I,lmax); eval_node3(Rz,fn,kiddie7,DIFF,I,lmax);
	}
}
void eval_node3(double (*fn)(double,double,double),node3 parent,location3 DIFF,
		double *I,double *J,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for(int i=0;i<8;i++) m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from  
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
	}else{
		eval_node3(fn,kiddie0,DIFF,I,J,lmax); eval_node3(fn,kiddie1,DIFF,I,J,lmax);
		eval_node3(fn,kiddie2,DIFF,I,J,lmax); eval_node3(fn,kiddie3,DIFF,I,J,lmax);
		eval_node3(fn,kiddie4,DIFF,I,J,lmax); eval_node3(fn,kiddie5,DIFF,I,J,lmax);
		eval_node3(fn,kiddie6,DIFF,I,J,lmax); eval_node3(fn,kiddie7,DIFF,I,J,lmax);
	}
}
void eval_node3(double *Rz,double (*fn)(double*,double,double,double),node3 parent,location3 DIFF,
		double *I,double *J,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for (int i=0;i<8;i++) m[i]=(*fn)(Rz,BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(Rz,BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from  
	f1=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(Rz,BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(Rz,BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(Rz,BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(Rz,BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
	}else{
		eval_node3(Rz,fn,kiddie0,DIFF,I,J,lmax); eval_node3(Rz,fn,kiddie1,DIFF,I,J,lmax);
		eval_node3(Rz,fn,kiddie2,DIFF,I,J,lmax); eval_node3(Rz,fn,kiddie3,DIFF,I,J,lmax);
		eval_node3(Rz,fn,kiddie4,DIFF,I,J,lmax); eval_node3(Rz,fn,kiddie5,DIFF,I,J,lmax);
		eval_node3(Rz,fn,kiddie6,DIFF,I,J,lmax); eval_node3(Rz,fn,kiddie7,DIFF,I,J,lmax);
	}
}
void eval_node3(double (*fn)(double,double,double),node3 parent,location3 DIFF,
		double *I,double *J,double *L,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,L1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for (int i=0;i<8;i++) m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from  
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
		L1=kiddie0.GetIysq()+kiddie1.GetIysq()+kiddie2.GetIysq()+kiddie3.GetIysq()+kiddie4.GetIysq()+
		   kiddie5.GetIysq()+kiddie6.GetIysq()+kiddie7.GetIysq();
		*L+=L1*dx*dy*dz;
	}else{
		eval_node3(fn,kiddie0,DIFF,I,J,L,lmax); eval_node3(fn,kiddie1,DIFF,I,J,L,lmax);
		eval_node3(fn,kiddie2,DIFF,I,J,L,lmax); eval_node3(fn,kiddie3,DIFF,I,J,L,lmax);
		eval_node3(fn,kiddie4,DIFF,I,J,L,lmax); eval_node3(fn,kiddie5,DIFF,I,J,L,lmax);
		eval_node3(fn,kiddie6,DIFF,I,J,L,lmax); eval_node3(fn,kiddie7,DIFF,I,J,L,lmax);
	}
}
void eval_node3(double *Rz,double (*fn)(double*,double,double,double),node3 parent,location3 DIFF,
		double *I,double *J,double *L,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,L1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for (int i=0;i<8;i++) m[i]=(*fn)(Rz,BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(Rz,BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from  
	f1=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(Rz,BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(Rz,BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(Rz,BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(Rz,BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
		L1=kiddie0.GetIysq()+kiddie1.GetIysq()+kiddie2.GetIysq()+kiddie3.GetIysq()+kiddie4.GetIysq()+
		   kiddie5.GetIysq()+kiddie6.GetIysq()+kiddie7.GetIysq();
		*L+=L1*dx*dy*dz;
	}else{
		eval_node3(Rz,fn,kiddie0,DIFF,I,J,L,lmax); eval_node3(Rz,fn,kiddie1,DIFF,I,J,L,lmax);
		eval_node3(Rz,fn,kiddie2,DIFF,I,J,L,lmax); eval_node3(Rz,fn,kiddie3,DIFF,I,J,L,lmax);
		eval_node3(Rz,fn,kiddie4,DIFF,I,J,L,lmax); eval_node3(Rz,fn,kiddie5,DIFF,I,J,L,lmax);
		eval_node3(Rz,fn,kiddie6,DIFF,I,J,L,lmax); eval_node3(Rz,fn,kiddie7,DIFF,I,J,L,lmax);
	}
}
void eval_node3(double (*fn)(double,double,double),node3 parent,location3 DIFF,
		double *I,double *J,double *K,double *L,double *M,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,K1,L1,M1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for (int i=0;i<8;i++) m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from  
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
		K1=kiddie0.GetIxsq()+kiddie1.GetIxsq()+kiddie2.GetIxsq()+kiddie3.GetIxsq()+kiddie4.GetIxsq()+
		   kiddie5.GetIxsq()+kiddie6.GetIxsq()+kiddie7.GetIxsq();
		*K+=K1*dx*dy*dz;
		L1=kiddie0.GetIysq()+kiddie1.GetIysq()+kiddie2.GetIysq()+kiddie3.GetIysq()+kiddie4.GetIysq()+
		   kiddie5.GetIysq()+kiddie6.GetIysq()+kiddie7.GetIysq();
		*L+=L1*dx*dy*dz;
		M1=kiddie0.GetIzsq()+kiddie1.GetIzsq()+kiddie2.GetIzsq()+kiddie3.GetIzsq()+kiddie4.GetIzsq()+
		   kiddie5.GetIzsq()+kiddie6.GetIzsq()+kiddie7.GetIzsq();
		*M+=M1*dx*dy*dz;
	}else{
		eval_node3(fn,kiddie0,DIFF,I,J,K,L,M,lmax); eval_node3(fn,kiddie1,DIFF,I,J,K,L,M,lmax);
		eval_node3(fn,kiddie2,DIFF,I,J,K,L,M,lmax); eval_node3(fn,kiddie3,DIFF,I,J,K,L,M,lmax);
		eval_node3(fn,kiddie4,DIFF,I,J,K,L,M,lmax); eval_node3(fn,kiddie5,DIFF,I,J,K,L,M,lmax);
		eval_node3(fn,kiddie6,DIFF,I,J,K,L,M,lmax); eval_node3(fn,kiddie7,DIFF,I,J,K,L,M,lmax);
	}
}
void eval_node3(double *Rz,double (*fn)(double*,double,double,double),node3 parent,location3 DIFF,
		double *I,double *J,double *K,double *L,double *M,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,K1,L1,M1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for (int i=0;i<8;i++) m[i]=(*fn)(Rz,BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(Rz,BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from  
	f1=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(Rz,BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(Rz,BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(Rz,BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(Rz,BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
		K1=kiddie0.GetIxsq()+kiddie1.GetIxsq()+kiddie2.GetIxsq()+kiddie3.GetIxsq()+kiddie4.GetIxsq()+
		   kiddie5.GetIxsq()+kiddie6.GetIxsq()+kiddie7.GetIxsq();
		*K+=K1*dx*dy*dz;
		L1=kiddie0.GetIysq()+kiddie1.GetIysq()+kiddie2.GetIysq()+kiddie3.GetIysq()+kiddie4.GetIysq()+
		   kiddie5.GetIysq()+kiddie6.GetIysq()+kiddie7.GetIysq();
		*L+=L1*dx*dy*dz;
		M1=kiddie0.GetIzsq()+kiddie1.GetIzsq()+kiddie2.GetIzsq()+kiddie3.GetIzsq()+kiddie4.GetIzsq()+
		   kiddie5.GetIzsq()+kiddie6.GetIzsq()+kiddie7.GetIzsq();
		*M+=M1*dx*dy*dz;
	}else{
		eval_node3(Rz,fn,kiddie0,DIFF,I,J,K,L,M,lmax); eval_node3(Rz,fn,kiddie1,DIFF,I,J,K,L,M,lmax);
		eval_node3(Rz,fn,kiddie2,DIFF,I,J,K,L,M,lmax); eval_node3(Rz,fn,kiddie3,DIFF,I,J,K,L,M,lmax);
		eval_node3(Rz,fn,kiddie4,DIFF,I,J,K,L,M,lmax); eval_node3(Rz,fn,kiddie5,DIFF,I,J,K,L,M,lmax);
		eval_node3(Rz,fn,kiddie6,DIFF,I,J,K,L,M,lmax); eval_node3(Rz,fn,kiddie7,DIFF,I,J,K,L,M,lmax);
	}
}
void eval_node3(double *Rz,double (*fn)(double*,double,double,double),node3 parent,location3 DIFF,
		double *I,double *J,double *K,double *L,double *M,double *N,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,K1,L1,M1,N1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for (int i=0;i<8;i++) m[i]=(*fn)(Rz,BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(Rz,BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from  
	f1=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(Rz,BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(Rz,BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(Rz,BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(Rz,BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb()+kiddie4.GetIyb()+
		   kiddie5.GetIyb()+kiddie6.GetIyb()+kiddie7.GetIyb();
		*J+=J1*dx*dy*dz;
		K1=kiddie0.GetIxsq()+kiddie1.GetIxsq()+kiddie2.GetIxsq()+kiddie3.GetIxsq()+kiddie4.GetIxsq()+
		   kiddie5.GetIxsq()+kiddie6.GetIxsq()+kiddie7.GetIxsq();
		*K+=K1*dx*dy*dz;
		L1=kiddie0.GetIysq()+kiddie1.GetIysq()+kiddie2.GetIysq()+kiddie3.GetIysq()+kiddie4.GetIysq()+
		   kiddie5.GetIysq()+kiddie6.GetIysq()+kiddie7.GetIysq();
		*L+=L1*dx*dy*dz;
		M1=kiddie0.GetIzsq()+kiddie1.GetIzsq()+kiddie2.GetIzsq()+kiddie3.GetIzsq()+kiddie4.GetIzsq()+
		   kiddie5.GetIzsq()+kiddie6.GetIzsq()+kiddie7.GetIzsq();
		*M+=M1*dx*dy*dz;
		N1=kiddie0.GetIxz()+kiddie1.GetIxz()+kiddie2.GetIxz()+kiddie3.GetIxz()+kiddie4.GetIxz()+
		   kiddie5.GetIxz()+kiddie6.GetIxz()+kiddie7.GetIxz();
		*N+=N1*dx*dy*dz;
	}else{
		eval_node3(Rz,fn,kiddie0,DIFF,I,J,K,L,M,N,lmax); eval_node3(Rz,fn,kiddie1,DIFF,I,J,K,L,M,N,lmax);
		eval_node3(Rz,fn,kiddie2,DIFF,I,J,K,L,M,N,lmax); eval_node3(Rz,fn,kiddie3,DIFF,I,J,K,L,M,N,lmax);
		eval_node3(Rz,fn,kiddie4,DIFF,I,J,K,L,M,N,lmax); eval_node3(Rz,fn,kiddie5,DIFF,I,J,K,L,M,N,lmax);
		eval_node3(Rz,fn,kiddie6,DIFF,I,J,K,L,M,N,lmax); eval_node3(Rz,fn,kiddie7,DIFF,I,J,K,L,M,N,lmax);
	}
}
void eval_node3(double *Rz,double (*fn)(double*,double,double,double),
		 double (*fL)(double),node3 parent,location3 DIFF,
		 double *I,double *J,double *K,double *L,double *M,double *N,int lmax){
	double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f14,f15,f16,f17,f18;
	double I0,I1,J1,K1,L1,M1,N1,dx,dy,dz,c0[9],m[8];
	location3 BottomLeft[8],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y; dz=.5*diff.z;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<8;i++){//bottom-left corners of children clockwise from origin
		BottomLeft[i].x=zero.x+fix3[i].x*dx;
		BottomLeft[i].y=zero.y+fix3[i].y*dy;
		BottomLeft[i].z=zero.z+fix3[i].z*dz;
	}
	for (int i=0;i<8;i++) m[i]=(*fn)(Rz,BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy,BottomLeft[i].z+.5*dz);
	f0=(*fn)(Rz,BottomLeft[1].x,BottomLeft[1].y,BottomLeft[1].z);//fn values clockwise from
	f1=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y+dy,BottomLeft[2].z);//origin
	f2=(*fn)(Rz,BottomLeft[2].x+dx,BottomLeft[2].y,BottomLeft[2].z);
	f3=(*fn)(Rz,BottomLeft[3].x,BottomLeft[3].y,BottomLeft[3].z);
	f4=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y,BottomLeft[2].z);
	f5=(*fn)(Rz,BottomLeft[4].x,BottomLeft[4].y,BottomLeft[4].z);//2nd layer starts
	f6=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z);
	f7=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y+dy,BottomLeft[5].z);
	f8=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z);
	f9=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y+dy,BottomLeft[6].z);
	f10=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z);
	f11=(*fn)(Rz,BottomLeft[7].x+dx,BottomLeft[7].y,BottomLeft[7].z);
	f12=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z);
	f14=(*fn)(Rz,BottomLeft[5].x,BottomLeft[5].y,BottomLeft[5].z+dz);//start top layer
	f15=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y+dy,BottomLeft[6].z+dz);
	f16=(*fn)(Rz,BottomLeft[6].x+dx,BottomLeft[6].y,BottomLeft[6].z+dz);
	f17=(*fn)(Rz,BottomLeft[7].x,BottomLeft[7].y,BottomLeft[7].z+dz);
	f18=(*fn)(Rz,BottomLeft[6].x,BottomLeft[6].y,BottomLeft[6].z+dz);
	node3 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,f4,f3,f5,f6,c0[8],f12,m[0]);
	node3 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,f4,f6,f7,f8,c0[8],m[1]);
	node3 kiddie2(BottomLeft[2],l1,DIFF,f4,f1,c0[2],f2,c0[8],f8,f9,f10,m[2]);
	node3 kiddie3(BottomLeft[3],l1,DIFF,f3,f4,f2,c0[3],f12,c0[8],f10,f11,m[3]);
	node3 kiddie4(BottomLeft[4],l1,DIFF,f5,f6,c0[8],f12,c0[4],f14,f18,f17,m[4]);
	node3 kiddie5(BottomLeft[5],l1,DIFF,f6,f7,f8,c0[8],f14,c0[5],f15,f18,m[5]);
	node3 kiddie6(BottomLeft[6],l1,DIFF,c0[8],f8,f9,f10,f18,f15,c0[6],f16,m[6]);
	node3 kiddie7(BottomLeft[7],l1,DIFF,f12,c0[8],f10,f11,f17,f18,f16,c0[7],m[7]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI()+kiddie4.GetI()+
	   kiddie5.GetI()+kiddie6.GetI()+kiddie7.GetI();
	if((l1>LMIN && fabs(8*I0-I1)<SMALL) || l1>lmax){
		*I+=I1*dx*dy*dz;
		J1=kiddie0.GetIfy(Rz,fL)+kiddie1.GetIfy(Rz,fL)+kiddie2.GetIfy(Rz,fL)
		   +kiddie3.GetIfy(Rz,fL)+kiddie4.GetIfy(Rz,fL)
		   +kiddie5.GetIfy(Rz,fL)+kiddie6.GetIfy(Rz,fL)+kiddie7.GetIfy(Rz,fL);
		*J+=J1*dx*dy*dz;
		K1=kiddie0.GetIxsq()+kiddie1.GetIxsq()+kiddie2.GetIxsq()+kiddie3.GetIxsq()+kiddie4.GetIxsq()+
		   kiddie5.GetIxsq()+kiddie6.GetIxsq()+kiddie7.GetIxsq();
		*K+=K1*dx*dy*dz;
		L1=kiddie0.GetIysq()+kiddie1.GetIysq()+kiddie2.GetIysq()+kiddie3.GetIysq()+kiddie4.GetIysq()+
		   kiddie5.GetIysq()+kiddie6.GetIysq()+kiddie7.GetIysq();
		*L+=L1*dx*dy*dz;
		M1=kiddie0.GetIzsq()+kiddie1.GetIzsq()+kiddie2.GetIzsq()+kiddie3.GetIzsq()+kiddie4.GetIzsq()+
		   kiddie5.GetIzsq()+kiddie6.GetIzsq()+kiddie7.GetIzsq();
		*M+=M1*dx*dy*dz;
		N1=kiddie0.GetIxz()+kiddie1.GetIxz()+kiddie2.GetIxz()+kiddie3.GetIxz()+kiddie4.GetIxz()+
		   kiddie5.GetIxz()+kiddie6.GetIxz()+kiddie7.GetIxz();
		*N+=N1*dx*dy*dz;
	}else{
		eval_node3(Rz,fn,fL,kiddie0,DIFF,I,J,K,L,M,N,lmax); eval_node3(Rz,fn,fL,kiddie1,DIFF,I,J,K,L,M,N,lmax);
		eval_node3(Rz,fn,fL,kiddie2,DIFF,I,J,K,L,M,N,lmax); eval_node3(Rz,fn,fL,kiddie3,DIFF,I,J,K,L,M,N,lmax);
		eval_node3(Rz,fn,fL,kiddie4,DIFF,I,J,K,L,M,N,lmax); eval_node3(Rz,fn,fL,kiddie5,DIFF,I,J,K,L,M,N,lmax);
		eval_node3(Rz,fn,fL,kiddie6,DIFF,I,J,K,L,M,N,lmax); eval_node3(Rz,fn,fL,kiddie7,DIFF,I,J,K,L,M,N,lmax);
	}
}
