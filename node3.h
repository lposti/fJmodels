#define LMIN 2

class location3{
	public:
		double x,y,z;
};

class intloc3{
	public:
		int x,y,z;
};

class node3{
	private:
		location3 loc;
		int level;
		double c_value[9];
		location3 diff;
	public:
		node3(void){};
		node3(double,double,double,int,location3);
		node3(location3,int,location3,double,double,double,double,double,double,double,double,double);
		~node3(void);
		int GetLevel(void){return level;};
		location3 GetLoc(void){return loc;};
		location3 Getdiff(void){return diff;};
		void GetVerts(location3 *);
		void GetCvals(double *);
		double GetSum(void);
		double GetVariance(double*);
		double GetI(void);
		double GetIyb(void);
		double GetIfy(double *,double (*)(double));
		double GetIxsq(void);
		double GetIysq(void);
		double GetIzsq(void);
		double GetIxz(void);
		void set_values(double (*)(double,double,double));
		void set_values(double*,double (*)(double*,double,double,double));
};
void eval_node3(double (*)(double,double,double),node3,location3,double*,int);
void eval_node3(double (*)(double,double,double),node3,location3,double*,double*,int);
void eval_node3(double (*)(double,double,double),node3,location3,double*,double*,double*,int);
void eval_node3(double (*)(double,double,double),node3,location3,double*,double*,double*,double*,double*,int);

void eval_node3(double*,double (*)(double*,double,double,double),node3,location3,double*,int);
void eval_node3(double*,double (*)(double*,double,double,double),node3,location3,double*,double*,int);
void eval_node3(double*,double (*)(double*,double,double,double),node3,location3,double*,double*,double*,int);
void eval_node3(double*,double (*)(double*,double,double,double),node3,location3,double*,double*,double*,double*,double*,int);
void eval_node3(double*,double (*)(double*,double,double,double),node3,location3,double*,double*,double*,double*,double*,double*,int);
void eval_node3(double*,double (*)(double*,double,double,double),double (*)(double),node3,location3,double*,double*,double*,double*,double*,double*,int);

