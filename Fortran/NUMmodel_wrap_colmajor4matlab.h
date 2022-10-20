void f_setupgeneralistsonly(const int n);

void f_setupgeneralistssimplepom(const int n, const int nPOM);

void f_setupgeneralistssimpleonly(const int n);

void f_setupgeneralistsonly_csp();

void f_setupdiatomsonly(const int n);

void f_setupdiatoms_simpleonly(const int n);

void f_setupgeneralistsdiatoms(const int n);

void f_setupgeneralistsdiatoms_simple(const int n);

void f_setupgeneralistssimplecopepod();

void f_setupgeneric(const int nCopepods, const double mAdult[]);

void f_setupnummodel(const int n, const int nCopepod, const int nPOM, const int nCopepods, const double mAdult[]);

void f_setupgeneric_csp(const int nCopepods, const double mAdult[]);

void f_sethtl(const double mHTL, const double mortHTL, const bool bQuadraticHTL, const bool bDecliningHTL);

void f_setmorthtl(const double mortHTL[]);

void f_calcderivatives(
		       const double u[],
		       const double L,
			   const double T,
		       const double dt,
		       double dudt[]);

/* void f_calcrates(
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
		       double g[]);*/

void f_simulatechemostateuler(
			      double u[],
			      const double L,
				  const double T,
				  const int nNutrients,
			      const double Ndeep[],
			      const double diff,
			      const double tEnd,
			      const double dt,
				  const bool bLosses);

void f_simulateeuler(
			      double u[],
			      const double L,
				  const double T, 
			      const double tEnd,
			      const double dt);

void f_getmass(
				double *m,
				double *mDelta);

void f_getsinking(double *velocity);

void f_getfunctions(
			double u[],
		    double *ProdGross,
		    double *ProdNet,
		    double *ProdHTL,
			double *ProdBact,
		    double *eHTL,
		    double *Bpico,
		    double *Bnano,
		    double *Bmicro);

void f_getbalance(
	const double u[],
	const double dudt[],
	double *Nbalance,
	double *Cbalance,
	double *Sibalance);
				  
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
		double *jPOM,
		double *mortpred,
		double *mortHTL,
		double *mort2,
		double *mort);
