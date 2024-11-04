%
% Get the HTL mortality and selection vectors
%
% In: 
%   p
% Out:
%   mortHTL - vector of the HTL mortality for all biomass groups
%   pHTL    - vector of the HTL selectivity for all biomass groups
%
function [mortHTL, pHTL] = getMortHTL(p)

mortHTL = zeros(1,p.n - p.idxB+1);
pHTL = mortHTL;
[mortHTL, pHTL] = calllib(loadNUMmodelLibrary(), 'f_getmorthtl', mortHTL, pHTL);

