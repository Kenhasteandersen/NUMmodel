%
% Returns the interaction matrix theta.
%
% In:
%  p - a structure with parameters from a setup
%
% Out:
%  theta - the size preference matrix of all size groups that are not
%          nutrients.
%
function theta = getTheta(p)

theta = zeros(p.n,p.n);
theta = calllib(loadNUMmodelLibrary, 'f_gettheta', theta);
theta = theta(p.idxB:end, p.idxB:end);

