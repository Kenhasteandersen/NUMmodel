function p = parametersGeneralists(n, mMax)
if (nargin==0)
    n = 25;
end
%
% Define parameters:
%
p.n = n;
[p.m, p.mLower, p.mDelta] = parametersCalcGrid(10^-8.5, mMax, p.n);

p.rhoCN = 5.68; % C:N mass ratio
p.epsilonL = 0.9; % Light uptake efficiency
p.epsilonF = 0.8; % Assimilation efficiency
p.cLeakage = 0.00015; % passive leakage of C and N
%
% Cell wall fraction of mass:
%
p.c = 0.0015; % the constant is increased a bit to limit the lower cell size
nu = p.c * p.m.^(-1/3);
%
% Clearance rates:
%
factor = (1e-6)^(1/3)/1.5;
p.AN = 0.00012; %0.00004 % 0.000162 % Mathilde.  (2.5e-3 l/d/cm) * (1e6 mug/g)^(-1/3) / 1.5 (g/cm); Andersen et al 2015
p.cN = 0.1;
p.AL = 0.000914; % if using Andys shading formula for non-diatoms
p.cL = 21; % if using Andys shading formula for non-diatoms
p.AF = 0.018;  %  Fits to TK data for protists
p.cF = 0.6; % Just a guess
%
% Calc rates as a function of m:
%
p.ANm = p.AN*p.m.^(1/3) ./ (1 + p.cN*p.m.^(-1/3));
p.ALm = p.AL*p.m.^(2/3) .* (1-exp(- p.cL*p.m.^(1/3) ));  % shading formula
p.AF = p.AF*p.m;
p.Jloss_passive_m = p.cLeakage * p.m.^(2/3); % in units of C
p.JFmax = p.cF*p.m.^(2/3);
%
% Prey encounter
%
p.beta = 500;
p.sigma = 1.3;
%
% Metabolism:
%
p.alphaJ = 1.5; % per day
p.Jmax = p.alphaJ * p.m .* (1-nu); % mugC/day
p.cR = 0.1;
p.Jresp = p.cR*p.alphaJ*p.m;
%
% Losses:
%
p.mort = 0*0.005*(p.Jmax./p.m) .* p.m.^(-1/4);
p.mort2 = 0.0002*double(p.n);

p.remin = 0.0; % fraction of mortality losses reminerilized to N and DOC
p.remin2 = 1; % fraction of virulysis remineralized to N and DOC
%
% Initial conditions:
%
p.N0 = 150;
p.DOC0 = 0;
p.B0 = ones(1,p.n);

end
