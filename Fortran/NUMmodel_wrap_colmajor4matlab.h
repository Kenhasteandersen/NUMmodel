void f_setupgeneralistsonly();

void f_setupgeneralistscopepod();

void f_setupgeneric(const int nCopepods, const double mAdult[]);

void f_calcderivatives(
		       const int nGrid,
		       const double u[],
		       const double L,
		       const double dt,
		       double dudt[]);

void f_calcrates(
		       const int nGrid,
		       const double u[],
		       const double L,
		       double jN[],
		       double jL[],
		       double jF[],
		       double jTot[],
		       double mortHTL[],
		       double g[]);

void f_simulatechemostateuler(
			      const int nGrid,
			      double u[],
			      const double L,
			      const double Ndeep,
			      const double diff,
			      const double tEnd,
			      const double dt); 

void f_simulateeuler(
			      const int nGrid,
			      double u[],
			      const double L,
			      const double tEnd,
			      const double dt); 
