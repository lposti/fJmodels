#include "math.h"
#include "uv_orb2.h"
#include "node1.h"


#define SMALL 1.e-7

double node1::GetI(void){
	return (c_value[0]+c_value[1]+4*c_value[2])/6.;// Simpson
}
void node1::GetCvals(double *vals){
	for(int i=0;i<3;i++) vals[i]=c_value[i];
}
void node1::GetVerts(double *verts){
	verts[0]=loc; verts[1]=loc+dx; verts[2]=loc+.5*dx;
}
void node1::set_values(double (*fn)(double)){
	double points[3];
	GetVerts(points);
	for(int i=0;i<3;i++) c_value[i]=(*fn)(points[i]);
}
void node1::set_values(uv_orb* orb,double (uv_orb::*fn)(double)){
	double points[3];
	GetVerts(points);
	for(int i=0;i<3;i++) c_value[i]=(orb->*fn)(points[i]);
}
node1::node1(double L,int lev,double DX,double value0,double value1,double value2){
	loc=L; level=lev; dx=DX/pow(2,level);
	c_value[0]=value0; c_value[1]=value1; c_value[2]=value2;
}
node1::node1(double lx,int lev,double DX){
	loc=lx; level=lev; dx=DX/pow(2,level);
}
void eval_node1(double (*fn)(double),node1 parent,double DX,double *I,int lmax){
	double I0,I1,m[2];
	double dx,c0[3];
	double BottomLeft[2],zero=parent.GetLoc();
	int l1=parent.GetLevel()+1;
	dx=.5*parent.Getdx();// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<2;i++){
		BottomLeft[i]=zero+i*dx;
		m[i]=(*fn)(BottomLeft[i]+.5*dx);
	}
	node1 kiddie0(BottomLeft[0],l1,DX,c0[0],c0[2],m[0]);
	node1 kiddie1(BottomLeft[1],l1,DX,c0[2],c0[1],m[1]);
	I1=kiddie0.GetI()+kiddie1.GetI();
	if(fabs(2*I0-I1)<SMALL || l1>lmax) *I+=I1*dx;
	else{
		eval_node1(fn,kiddie0,DX,I,lmax);
		eval_node1(fn,kiddie1,DX,I,lmax);
	}
}
void eval_node1(uv_orb *orb,double (uv_orb::*fn)(double),node1 parent,double DX,double *I,int lmax){
	double I0,I1,m[2];
	double dx,c0[3];
	double BottomLeft[2],zero=parent.GetLoc();
	int l1=parent.GetLevel()+1;
	dx=.5*parent.Getdx();// dx now half side length of parent
	I0=parent.GetI();
	parent.GetCvals(c0);
	for(int i=0;i<2;i++){
		BottomLeft[i]=zero+i*dx;
		m[i]=(orb->*fn)(BottomLeft[i]+.5*dx);
	}
	node1 kiddie0(BottomLeft[0],l1,DX,c0[0],c0[2],m[0]);
	node1 kiddie1(BottomLeft[1],l1,DX,c0[2],c0[1],m[1]);
	I1=kiddie0.GetI()+kiddie1.GetI();
	if(fabs(2*I0-I1)<SMALL || l1>lmax) *I+=I1*dx;
	else{
		eval_node1(orb,fn,kiddie0,DX,I,lmax);
		eval_node1(orb,fn,kiddie1,DX,I,lmax);
	}
}
node1::~node1(void){}
