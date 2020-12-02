void f_setparameters(
    const int n, 
    const double m[],
    const double rhoCN, 
    const double epsilonL,
    const double epsilonF,
    const double AN[],
    const double AL[],
    const double AF[],
    const double Jmax[], 
    const double JFmax[], 
    const double Jresp[],
    const double Jloss_passive[],
    const double* theta,
    const double mort[],
    const double mort2,
    const double mortHTL[],
    const double remin,
    const double remin2,
    const double cLeakage
);

 void f_calcrates(double T, double L, int n, double u[], double gammaN, double gammaDOC, double *dudt);
