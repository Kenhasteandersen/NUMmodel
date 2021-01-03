void f_setupgeneralistsonly();

void f_setupgeneralistscopepod();

void f_setupgeneric(const int nCopepods, double mAdult[]);

void f_calcderivatives(
		       const int nGrid,
		       const double u[],
		       const double L,
		       const double dt,
		       double dudt[]);

void f_simulatechemostateuler(
			      const int nGrid,
			      const int iEnd,
			      double usave[],
			      const double u0[],
			      const doubleL,
			      const double diff,
			      const double tEnd,
			      const double dt); 
