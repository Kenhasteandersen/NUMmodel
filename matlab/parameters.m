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
%==========================================
% Set up the specific systems to simulate
% =========================================

n = 10; % No of bins in each group
%
% Universal parameters:
%
p.nGroups = 2;
p.typeGroups = [1,2]; % 1=generalists, 2=copepods
%
% Parameters for unicellular generalists
%
p.pGeneralists = parametersGeneralists(n);
p.ixStart = 3;
p.ixEnd = p.ixStart+n-1;
%
% Parameters for copepods:
%
p.pCopepods = parametersCopepods(mCopepods,n);
p.ixStart(2) = p.ixEnd(1)+1;
p.ixEnd(2) = p.ixStart(2)+n-1;
% =============================================
% Set up derived parameters:
% =============================================
% Mass grid:
p.m = [NaN, NaN, p.pGeneralists.m, p.pCopepods.m];
p.mLower = [NaN, NaN, p.pGeneralists.mLower, p.pCopepods.mLower];
p.mDelta = [NaN, NaN, p.pGeneralists.mDelta, p.pCopepods.mDelta];
% Feeding parameters:
p.AF = [NaN, NaN, p.pGeneralists.AFm, p.pCopepods.AF];
p.JFmax = [NaN, NaN, p.pGeneralists.JFmaxm, p.pCopepods.JFmax];
p.epsilonF(p.ixStart(1):p.ixEnd(1)) = p.pGeneralists.epsilonF;
p.epsilonF(p.ixStart(2):p.ixEnd(2)) = p.pCopepods.epsilonF;
p.Jresp(p.ixStart(1):p.ixEnd(1)) = p.pGeneralists.Jrespm;
p.Jresp(p.ixStart(2):p.ixEnd(2)) = p.pCopepods.Jresp;
%
% Calculate interaction matrix:
%
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
p.mortHTLm = 0.1*(1./(1+(p.m./p.mHTL).^(-2)));
% ================================================================
%  Initial conditions (also used for deep conditions of chemostat
% ================================================================
p.u0 = 10*p.m./p.m;
p.u0(1) = 150;
p.u0(2) = 0;

