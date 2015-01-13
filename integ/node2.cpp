#include <math.h>
#include "node2.h"

#define SMALL 1.e-7
intloc2 fix2[4]={{0,0},{0,1},{1,1},{1,0}};

double node2::GetSum(void){
	double sum=c_value[0];
	for(int i=1;i<4;i++) sum+=c_value[i];
	return sum;
}
double node2::GetVariance(double *sum1){
	double sum2=pow(c_value[0],2); *sum1=c_value[0];
	for(int i=1;i<4;i++){
		*sum1+=c_value[i]; sum2+=pow(c_value[i],2);
	}
	return .25*(sum2-.25*pow(*sum1,2));
}
double node2::GetI(void){
	double sum=c_value[0];
	for(int i=1;i<4;i++) sum+=c_value[i];
	return .5*(c_value[4]+.25*sum);
}
double node2::GetIyb(void){//for evaluation of rho*ybar 
	double sum=0;
	for(int i=0;i<4;i++) sum+=c_value[i]*(loc.y+fix2[i].y*diff.y);
	return .5*(c_value[4]*(loc.y+0.5*diff.y)+.25*sum);
}
void node2::GetCvals(double *vals){
	for(int i=0;i<5;i++) vals[i]=c_value[i];
}
void node2::GetVerts(location2 *verts){
	for(int i=0;i<4;i++){
		verts[i].x=loc.x+fix2[i].x*diff.x;
		verts[i].y=loc.y+fix2[i].y*diff.y;
	}
	verts[4].x=loc.x+.5*diff.x; verts[4].y=loc.y+.5*diff.y;
}
void node2::set_values(double (*fn)(double,double)){
	location2 points[5];
	GetVerts(points);
	for(int i=0;i<5;i++) c_value[i]=(*fn)(points[i].x,points[i].y);
}
void node2::set_values(double *Rz,double (*fn)(double*,double,double)){
	location2 points[5];
	GetVerts(points);
	for(int i=0;i<5;i++) c_value[i]=(*fn)(Rz,points[i].x,points[i].y);
}
node2::node2(location2 L,int lev,location2 DIFF,double value0,double value1,double value2,
	     double value3,double value4){
	loc=L; level=lev; double k=pow(2,level);
	diff.x=DIFF.x/k; diff.y=DIFF.y/k;
	c_value[0]=value0; c_value[1]=value1; c_value[2]=value2;
	c_value[3]=value3; c_value[4]=value4;
}
node2::node2(double lx,double ly,int lev,location2 DIFF){
	loc.x=lx; loc.y=ly; level=lev; double k=pow(2,level);
	diff.x=DIFF.x/k; diff.y=DIFF.y/k;
}
node2::~node2(void){}
void eval_node2(double (*fn)(double,double),node2 parent,location2 DIFF,double *I,int lmax){
	double f0,f1,f2,f3,I0,I1,dx,dy,c0[5],m[4];
	location2 BottomLeft[4],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<4;i++){//bottom-left corners & central values of children
		BottomLeft[i].x=zero.x+fix2[i].x*dx;// clockwise from origin
		BottomLeft[i].y=zero.y+fix2[i].y*dy;
		m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy);
	}
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y);   //fn values clockwise from  
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y);
	node2 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,c0[4],f3,m[0]);
	node2 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,c0[4],m[1]);
	node2 kiddie2(BottomLeft[2],l1,DIFF,c0[4],f1,c0[2],f2,m[2]);
	node2 kiddie3(BottomLeft[3],l1,DIFF,f3,c0[4],f2,c0[3],m[3]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI();
//	if(l1>3 && fabs(4*I1-I0)<SMALL) I+=I1*dx*dy;
//	if(l1>7) printf("%d %g (%g,%g) ",l1,4*I1-I0,zero.x,zero.y);
	if((l1>LMIN && fabs(4*I0-I1)<.25*SMALL) || l1>lmax) *I+=I1*dx*dy;
	else{
		eval_node2(fn,kiddie0,DIFF,I,lmax); eval_node2(fn,kiddie1,DIFF,I,lmax);
		eval_node2(fn,kiddie2,DIFF,I,lmax); eval_node2(fn,kiddie3,DIFF,I,lmax);
	}
}
void eval_node2(double *Rz,double (*fn)(double*,double,double),node2 parent,location2 DIFF,double *I,int lmax){
	double f0,f1,f2,f3,I0,I1,dx,dy,c0[5],m[4];
	location2 BottomLeft[4],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<4;i++){//bottom-left corners & central values of children
		BottomLeft[i].x=zero.x+fix2[i].x*dx;// clockwise from origin
		BottomLeft[i].y=zero.y+fix2[i].y*dy;
		m[i]=(*fn)(Rz,BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy);
	}
	f0=(*fn)(Rz,BottomLeft[1].x,BottomLeft[1].y);   //fn values clockwise from  
	f1=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y+dy);//origin
	f2=(*fn)(Rz,BottomLeft[2].x+dx,BottomLeft[2].y);
	f3=(*fn)(Rz,BottomLeft[3].x,BottomLeft[3].y);
	node2 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,c0[4],f3,m[0]);
	node2 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,c0[4],m[1]);
	node2 kiddie2(BottomLeft[2],l1,DIFF,c0[4],f1,c0[2],f2,m[2]);
	node2 kiddie3(BottomLeft[3],l1,DIFF,f3,c0[4],f2,c0[3],m[3]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI();
//	if(l1>3 && fabs(4*I1-I0)<SMALL) I+=I1*dx*dy;
//	if(l1>7) printf("%d %g (%g,%g) ",l1,4*I1-I0,zero.x,zero.y);
	if((l1>LMIN && fabs(4*I0-I1)<.25*SMALL) || l1>lmax) *I+=I1*dx*dy;
	else{
		eval_node2(Rz,fn,kiddie0,DIFF,I,lmax); eval_node2(Rz,fn,kiddie1,DIFF,I,lmax);
		eval_node2(Rz,fn,kiddie2,DIFF,I,lmax); eval_node2(Rz,fn,kiddie3,DIFF,I,lmax);
	}
}
void eval_node2(double (*fn)(double,double),node2 parent,location2 DIFF,
		double *I,double *J,int lmax){
	double f0,f1,f2,f3,I0,I1,J1,dx,dy,c0[5],m[4];
	location2 BottomLeft[4],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<4;i++){//bottom-left corners & central values of children
		BottomLeft[i].x=zero.x+fix2[i].x*dx;// clockwise from origin
		BottomLeft[i].y=zero.y+fix2[i].y*dy;
		m[i]=(*fn)(BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy);
	}
	f0=(*fn)(BottomLeft[1].x,BottomLeft[1].y);   //fn values clockwise from  
	f1=(*fn)(BottomLeft[2].x,BottomLeft[2].y+dy);//origin
	f2=(*fn)(BottomLeft[2].x+dx,BottomLeft[2].y);
	f3=(*fn)(BottomLeft[3].x,BottomLeft[3].y);
	node2 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,c0[4],f3,m[0]);
	node2 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,c0[4],m[1]);
	node2 kiddie2(BottomLeft[2],l1,DIFF,c0[4],f1,c0[2],f2,m[2]);
	node2 kiddie3(BottomLeft[3],l1,DIFF,f3,c0[4],f2,c0[3],m[3]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI();
	if((l1>LMIN && fabs(4*I0-I1)<.25*SMALL) || l1>lmax){
		*I+=I1*dx*dy;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb();
		*J+=J1*dx*dy;
	}else{
		eval_node2(fn,kiddie0,DIFF,I,J,lmax); eval_node2(fn,kiddie1,DIFF,I,J,lmax);
		eval_node2(fn,kiddie2,DIFF,I,J,lmax); eval_node2(fn,kiddie3,DIFF,I,J,lmax);
	}
}
void eval_node2(double *Rz,double (*fn)(double*,double,double),node2 parent,location2 DIFF,
		double *I,double *J,int lmax){
	double f0,f1,f2,f3,I0,I1,J1,dx,dy,c0[5],m[4];
	location2 BottomLeft[4],zero=parent.GetLoc(),diff=parent.Getdiff();
	int l1=parent.GetLevel()+1;
	dx=.5*diff.x; dy=.5*diff.y;// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<4;i++){//bottom-left corners & central values of children
		BottomLeft[i].x=zero.x+fix2[i].x*dx;// clockwise from origin
		BottomLeft[i].y=zero.y+fix2[i].y*dy;
		m[i]=(*fn)(Rz,BottomLeft[i].x+.5*dx,BottomLeft[i].y+.5*dy);
	}
	f0=(*fn)(Rz,BottomLeft[1].x,BottomLeft[1].y);   //fn values clockwise from  
	f1=(*fn)(Rz,BottomLeft[2].x,BottomLeft[2].y+dy);//origin
	f2=(*fn)(Rz,BottomLeft[2].x+dx,BottomLeft[2].y);
	f3=(*fn)(Rz,BottomLeft[3].x,BottomLeft[3].y);
	node2 kiddie0(BottomLeft[0],l1,DIFF,c0[0],f0,c0[4],f3,m[0]);
	node2 kiddie1(BottomLeft[1],l1,DIFF,f0,c0[1],f1,c0[4],m[1]);
	node2 kiddie2(BottomLeft[2],l1,DIFF,c0[4],f1,c0[2],f2,m[2]);
	node2 kiddie3(BottomLeft[3],l1,DIFF,f3,c0[4],f2,c0[3],m[3]);
	I1=kiddie0.GetI()+kiddie1.GetI()+kiddie2.GetI()+kiddie3.GetI();
	if((l1>LMIN && fabs(4*I0-I1)<.25*SMALL) || l1>lmax){
		*I+=I1*dx*dy;
		J1=kiddie0.GetIyb()+kiddie1.GetIyb()+kiddie2.GetIyb()+kiddie3.GetIyb();
		*J+=J1*dx*dy;
	}else{
		eval_node2(Rz,fn,kiddie0,DIFF,I,J,lmax); eval_node2(Rz,fn,kiddie1,DIFF,I,J,lmax);
		eval_node2(Rz,fn,kiddie2,DIFF,I,J,lmax); eval_node2(Rz,fn,kiddie3,DIFF,I,J,lmax);
	}
}
