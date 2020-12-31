function p = parametersInit()
%
% Constants used for indexing etc.:
%
p.ixN = 1;
p.ixDOC = 2;
p.typeGeneralists = 1;
p.typeCopepods = 10;
% =============================================
% Set up dummy parameters for nutrient fields:
% =============================================
% Mass grid:
p.m = [NaN, NaN];
p.mLower = [NaN, NaN];
p.mDelta = [NaN, NaN];
% Feeding parameters:
p.AF = [NaN, NaN];
p.JFmax = [NaN, NaN];
p.epsilonF = [NaN, NaN];
p.Jresp = [NaN, NaN];

p.nGroups = 0;