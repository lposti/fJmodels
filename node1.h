class node1{
	private:
		double loc;
		int level;
		double c_value[3];
		double dx;
	public:
		node1(void){};
		node1(double,int,double);
		node1(double,int,double,double,double,double);
		~node1(void);
		int GetLevel(void){return level;};
		double GetLoc(void){return loc;};
		double Getdx(void){return dx;};
		void GetVerts(double *);
		void GetCvals(double *);
		double GetI(void);
		void set_values(double (*)(double));
		void set_values(uv_orb *orb,double (uv_orb::*)(double));
};
void eval_node1(double (*)(double),node1,double,double*,int);
void eval_node1(uv_orb*,double (uv_orb::*)(double),node1,double,double*,int);
