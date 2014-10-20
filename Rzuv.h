#include <math.h>

class pt{
	private:
		double R2;
		double z2;
		double r2;
	public:
		pt(double,double,double);
		pt(double,double*,double*);
		double Delta;
		double Delta2;
		double R;
		double z;
		double u;
		double v;
		double shu;
		double chu;
		double shu2;
		double chu2;
		double sv;
		double cv;
		double sv2;
		double cv2;
		double pR;
		double pz;
		double p2;
		double pu;
		double pv;
};
