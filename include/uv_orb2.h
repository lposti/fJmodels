class uv_orb{
	private:
		double Rc1;
		double Ju1;
		double Jv1;
		double I3u;
		double I3v;
		double JuE;
		double JuI;
		double JuL;
		double JvE;
		double JvI;
		double JvL;
		double psivmax;
		double Deltau;
		double Lzsq;
		double shu2;
		double sh1sq;
		double Phiu1;
		double ubar;
		double sv;
		double cv;
		double sv2;
		void GetTurnu(double);
		void GetJu(void);
		double GetTurnv(double);
		void GetJv(void);
		void Getu0(double);
		void GetRc(double);
		void GetFreqs(void);
		void infinity();
//		double Xbrent(uv_orb*,double (uv_orb::*)(double),double,double,double,double,double,int);
	public:
		uv_orb(void){};
		uv_orb(double,double,double,double*,double*);
		~uv_orb(void){};
		void reset(double,double,double,double*,double*);
		double uinner;
		double uouter;
		double umid;//minimum of effective pot
		double vturn;
		//double RcE(void);
		double Rc(void);
		double Ju(void);
		double Jv(void);
		void GetThetas(double*,double*,double*,double*,double*);
		double Delta;
		double u;
		double v;
		double pv;
		double Iu;
		double Iv;
		double E;
		double Lz;
		double Omegau;
		double Omegav;
		double Omegaphi;
		double Phiu(double);
		double dU(double);
		double dPhiu(double,double*);
		double uturnfn(double);
		double duturnfn(double);
		double Rcfn(double);
		double Juint(double);
		double dpudE(double);
		double dpudI3(double);
		double dpudLz(double);
		double Phiv(double);
		double dV(double);
		double vturnfn(double);
		double Jvint(double);
		double dpvdE(double);
		double dpvdI3(double);
		double dpvdLz(double);
};
