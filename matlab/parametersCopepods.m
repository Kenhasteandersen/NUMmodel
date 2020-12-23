function p = parametersCopepods(mAdult,n)
%
% Grid:
%
%       m(1)                           mAdult=m(end)
% |--------------|--  ...  -----|--------------|
% mLower(1)                    
% =mOffspring
%
% Offspring size is the lower size of the first grid cell.
% Adult size is the middle size of the last grid cell.
%
% Parameters for active copepods from Serra-Pompei et al (2020):
% (note: rates are absolute and not relative)
%
p.AdultOffspring = 100;
p.beta = 10000;
p.sigma = 1.5;
p.v = 0.011;
p.q = 0.75;
p.h = 1.37;
p.hExponent = 0.75;
p.Kappa = 0.16;
p.p = 0.75;
p.epsilonF = 0.67;
p.epsilonR = 0.25;
%
% Make grid:
%
p.n = n;
lnDelta = (log(mAdult)-log(mAdult/p.AdultOffspring)) / (n-0.5);
mMin = exp(log(mAdult/p.AdultOffspring)+0.5*lnDelta);
[p.m, p.mLower, p.mDelta, p.mZ] = parametersCalcGrid(mMin, mAdult, p.n);
%
% Derived parameters:
%
p.AF = p.v*p.m.^p.q; % Clearance rate
p.JFmax = p.h*p.m.^p.hExponent;
p.Jresp = p.Kappa*p.m.^p.hExponent;

