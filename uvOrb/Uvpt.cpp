/*
 * Uvpt.cpp
 *
 *  Created on: Feb 22, 2015
 *      Author: L. Posti
 */

#include <math.h>
#include "Uvpt.h"

pt::pt(const double Delta0,const double *x,const double *p){
	Delta=Delta0; Delta2=Delta*Delta;
	R=x[0]; z=x[1]; R2=R*R; z2=z*z; r2=R2+z2;
	shu2 =.5*(r2-Delta2+sqrt(pow(Delta2-r2,2)+4*Delta2*R2))/Delta2;
	shu=sqrt(shu2); u=asinh(shu);
	chu2=1+shu2; chu=sqrt(chu2);
	cv2=.5*(r2+Delta2-sqrt(pow(Delta2+r2,2)-4*Delta2*z2))/Delta2;
	cv= z<0? -sqrt(cv2) : sqrt(cv2);
	v=acos(cv); sv=sin(v);
	sv2=1-cv2;
	pR=p[0]; pz=p[1]; p2=pR*pR+pz*pz;
	pu=Delta*(pR*chu*sv+pz*shu*cv);
	pv=Delta*(pR*shu*cv-pz*chu*sv);
}
pt::pt(const double Delta0,const double ui,const double vi){
	Delta=Delta0; u=ui; v=vi;
	R=Delta*sinh(u)*sin(v);
	z=Delta*cosh(u)*cos(v);

	Delta2=Delta*Delta;
	shu=chu=shu2=chu2=0;
	sv=cv=sv2=cv2=0;
	pR=pz=p2=pu=pv=0;
	R2=z2=r2=0;
}

