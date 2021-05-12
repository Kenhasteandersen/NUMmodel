%
% Get the mass grid from the Fortran library and return it in the
% p-structure
%
function p = getMass(p)
%
% Get mass and cell width (Delta-mass) from the library:
%
p.m = zeros(1,p.n);
p.mDelta = p.m;
[p.m, p.mDelta] = calllib(loadNUMmodelLibrary(), 'f_getmass', p.m, p.mDelta);
%
% Calculate the upper and lower mass of a cell:
%
p.mLower = (sqrt(4*p.m.^2+p.mDelta.^2)-p.mDelta)/2;
p.mUpper = p.mLower+p.mDelta;

