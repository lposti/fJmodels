#define LMIN 2

class location2{
	public:
		double x,y;
};

class intloc2{
	public:
		int x,y;
};

class node2{
	private:
		location2 loc;
		int level;
		double c_value[5];
		location2 diff;
	public:
		node2(void){};
		node2(double,double,int,location2);
		node2(location2,int,location2,double,double,double,double,double);
		~node2(void);
		int GetLevel(void){return level;};
		location2 GetLoc(void){return loc;};
		location2 Getdiff(void){return diff;};
		void GetVerts(location2 *);
		void GetCvals(double *);
		double GetSum(void);
		double GetVariance(double*);
		double GetI(void);
		double GetIyb();
		void set_values(double (*)(double,double));
		void set_values(double*,double (*)(double*,double,double));
};
void eval_node2(double (*)(double,double),node2,location2,double*,int);
void eval_node2(double (*)(double,double),node2,location2,double*,double*,int);

void eval_node2(double*,double (*)(double*,double,double),node2,location2,double*,int);
void eval_node2(double*,double (*)(double*,double,double),node2,location2,double*,double*,int);
