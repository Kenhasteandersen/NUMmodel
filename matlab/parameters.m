function p = parameters(mCopepods)
%
% Clarification of the grid:
%
%   |           i          |
% --|----------------------|----
% mLower_i     m_i         
%
% m_i represents the geometric center of a cell
%
iGeneralists = 1;
iCopepods = 10;

n = 10; % No of bins in each group
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
p.Jresp = [NaN, NaN];
% =========================================
% Set up the specific systems to simulate
% =========================================
p.nGroups = 0;
%
% Parameters for unicellular generalists
%
p = parametersAddgroup(iGeneralists, p, n);
%
% Parameters for copepods:
%
for i = 1:length(mCopepods)
    p = parametersAddgroup(iCopepods, p, n, mCopepods(i));
end
% =========================================
% Calculate interaction matrix:
% =========================================
beta(p.ixStart(1):p.ixEnd(1)) = p.pGeneralists.beta;
sigma(p.ixStart(1):p.ixEnd(1)) = p.pGeneralists.sigma;

beta(p.ixStart(2):p.ixEnd(2)) = p.pCopepods.beta;
sigma(p.ixStart(2):p.ixEnd(2)) = p.pCopepods.sigma;

ix = 3:length(p.m);
p.theta = parametersCalcTheta(p.m(ix),p.mLower(ix),p.mDelta(ix),beta(ix),sigma(ix));
%
% Set higher trophic level (HTL) mortality:
%
betaHTL = 500; % The predator-prey size ratio of HTLs
p.mHTL = max(p.m)/betaHTL^1.5; % Bins affected by HTL mortality
p.mortHTLm = 0.01*(1./(1+(p.m./p.mHTL).^(-2)));
% ================================================================
%  Initial conditions (also used for deep conditions of chemostat
% ================================================================
p.u0 = 10*p.m./p.m; % Biomasses
p.u0(1) = 150; % Nutrients
p.u0(2) = 0; % DOC

