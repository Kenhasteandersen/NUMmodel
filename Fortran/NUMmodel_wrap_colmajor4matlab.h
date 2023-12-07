void f_setupgeneralistsonly(const int n, bool *Clost, char *errorstr);

void f_setupgeneralistssimplepom(const int n, const int nPOM, bool *Clost, char *errorstr);

void f_setupgeneralistspom(const int n, const int nPOM, bool *Clost, char *errorstr);

void f_setupgeneralistssimpleonly(const int n, bool *Clost, char *errorstr);

void f_setupdiatomsonly(const int n, bool *Clost, char *errorstr);

void f_setupdiatoms_simpleonly(const int n, bool *Clost, char *errorstr);

void f_setupgeneralistsdiatoms(const int n, bool *Clost, char *errorstr);

void f_setupgeneralistsdiatoms_simple(const int n, bool *Clost, char *errorstr);

void f_setupgeneralistssimplecopepod(bool *Clost, char *errorstr);

void f_setupgeneric(const int nCopepods, const double mAdult[], bool *Clost, char *errorstr);

void f_setupnummodel(const int n, const int nCopepod, const int nPOM, 
					const int nCopepodsPassive, const double mAdultPassive[], 
					const int nCopepodsActive, const double mAdultActive[], 
					bool *Clost, char *errorstr);

void f_setupnummodelsimple(const int n, const int nCopepod, const int nPOM, const int nCopepods, 
                                        const double mAdult[], bool *Clost, char *errorstr);

void f_setupgendiatcope(const int n,const int nCopepod, const int nPOM, const int nCopepods, 
                                        const double mAdult[], bool *Clost, char *errorstr);

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

void f_simulateeulerfunctions(
			      double u[],
			      const double L,
				  const double T, 
			      const double tEnd,
			      const double dt,
				  double *ProdGross,
				  double *ProdNet,
				  double *ProdHTL,
				  double *ProdBact,
				  double *eHTL,
				  double *Bpico,
				  double *Bnano,
				  double *Bmicro);

void f_getmass(
				double *m,
				double *mDelta);

void f_getsinking(double *velocity);

void f_setsinking(double *velocity);

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
	double *Cbalance,
	double *Nbalance,
	double *Sibalance);

void f_getlost(
	const double u[],
	double *Clost,
	double *Nlost,
	double *SiLost);
				  
void f_getrates(
		double *jN,
		double *jDOC,
		double *jL,
		double *jSi,
		double *jF,
		double *jFreal,
		double *f,
		double *jTot,
		double *jMax,
		double *jFmaxx,
		double *jR,
		double *jResptot,
		double *jLossPassive, 
		double *jNloss,
		double *jLreal, 
		double *jPOM,
		double *mortpred,
		double *mortHTL,
		double *mort2,
		double *mort);

void f_gettheta(
		double *thetaMatrix
);
