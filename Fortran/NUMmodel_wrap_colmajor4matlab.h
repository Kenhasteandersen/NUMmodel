void f_setupgeneralistsonly();

void f_setupgeneralistscopepod();

void f_setupgeneric(const int nCopepods, double mAdult[]);

void f_calcderivatives(
		       const int nGrid,
		       const double u[],
		       const double L,
		       const double dt,
		       double dudt[]);
