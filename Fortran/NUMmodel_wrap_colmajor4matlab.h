void f_setupgeneralistsonly(const int n);

void f_setupgeneralistsonly_csp();

void f_setupdiatomsonly(const int n);

void f_setupdiatoms_simpleonly(const int n);

void f_setupgeneralistsdiatoms(const int n);

void f_setupgeneralistsdiatoms_simple(const int n);

void f_setupgeneralistscopepod();

void f_setupgeneric(const int nCopepods, const double mAdult[]);

void f_setupgeneric_csp(const int nCopepods, const double mAdult[]);

void f_parametershtl(const double mortHTL, const bool bQuadraticHTL);

void f_calcderivatives(
		       const int nGrid,
		       const double u[],
		       const double L,
			   const double T,
		       const double dt,
		       double dudt[]);

void f_calcrates(
		       const int nGrid,
		       const double u[],
		       const double L,
			   const double T,
		       double jN[],
		       double jL[],
		       double jF[],
		       double jTot[],
		       double mortHTL[],
		       double mortpred[],
		       double g[]);

void f_simulatechemostateuler(
			      double u[],
			      const double L,
				  const double T,
				  const int nNutrients,
			      const double Ndeep[],
			      const double diff,
			      const double tEnd,
			      const double dt);

void f_simulateeuler(
			      const int nGrid,
			      double u[],
			      const double L,
				  const double T, 
			      const double tEnd,
			      const double dt);

void f_getmass(
				double *m,
				double *mDelta);

void f_getfunctions(
		    double *ProdGross,
		    double *ProdNet,
		    double *ProdHTL,
		    double *eHTL,
		    double *Bpico,
		    double *Bnano,
		    double *Bmicro);

void f_getbalance(double *Nbalance,
				  double *Cbalance);
				  
void f_getrates(
		double *jN,
		double *jDOC,
		double *jL,
		double *jSi,
		double *jF,
		double *jFreal,
		double *jTot,
		double *jMax,
		double *jFmaxx,
		double *jR,
		double *jLossPassive, 
		double *jNloss,
		double *jLreal, 
		double *mortpred,
		double *mortHTL,
		double *mort2,
		double *mort);
