function [Brange,f] = calcBiomassRange(B,m, m1, m2)
m2 = min(max(m),m2);
%
% Find lower and upper cell boundaries:
%
Delta = m(2)/m(1);
mLower = m/sqrt(Delta);
mUpper = m*sqrt(Delta);
%
% Find affected size ranges:
%
ix = find(mUpper>=m1 & mLower<=m2);
%
% Find fractions to take from boundary size groups:
%
fLower = 1-(log(m1)-log(mLower(ix(1)))) / (log(mUpper(ix(1)))-log(mLower(ix(1))));
fUpper = (log(m2)-log(mLower(ix(end)))) / (log(mUpper(ix(end)))-log(mLower(ix(end))));
%
% Assemble fractions:
%
f = zeros(1,length(m));
f(ix) = 1;
f(ix(1)) = fLower;
f(ix(end)) = fUpper;
%
% Assemble biomass
%
Brange = sum(B.*f); %fLower*B(ix(1)) + sum(B(ix(2:end-1))) + fUpper*B(ix(end));
%Brange = fLower*B(ix(1)) + sum(B(ix(2:end-1))) + fUpper*B(ix(end));
