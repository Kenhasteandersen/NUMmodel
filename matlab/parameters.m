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

p = parametersInit();
n = 10; % No of bins in each group
% =========================================
% Set up the specific systems to simulate
% =========================================

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
ix = 3:length(p.m);
p.theta = parametersCalcTheta(p.m(ix),p.mLower(ix),p.mDelta(ix),p.beta(ix),p.sigma(ix));
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

