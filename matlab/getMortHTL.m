%
% Get the HTL mortality and selection vectors
%
% In: 
%   p
% Out:
%   mortHTL - vector of the HTL mortality for all biomass groups. Must be
%       multiplied by B if bQuadratic = true.
%   bQuadratic - whether the mortality should be multiplied by B or not
%       (boolean)
%
function [mortHTL, bQuadratic] = getMortHTL(p)

mortHTL = zeros(1,p.n - p.idxB+1);
bQuadratic = false;
[mortHTL, bQuadratic] = calllib(loadNUMmodelLibrary(), 'f_getmorthtl', mortHTL, logical(bQuadratic));

