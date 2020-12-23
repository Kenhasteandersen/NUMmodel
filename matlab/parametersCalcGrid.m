%
% Make a grid. 
% Input: mMin and mMax are _center_ masses of the grid
%
function [m, mLower, mDelta, z] = parametersCalcGrid(mMin, mMax, n)

x = linspace(log(mMin), log(mMax), n);
deltax = x(2)-x(1);
m = exp(x);
mLower = exp(x-0.5*deltax);
mDelta =  exp(x+0.5*deltax)-mLower;
z = mLower./(mLower+mDelta);
