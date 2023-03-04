%
% Return the biomass within the size range m1:m2.
% Also output the fraction of biomass from each size
% group that falls into the range.
%
function [Brange,f] = calcBiomassRange(B,m, m1, m2)
B = reshape(B,1,length(B));

%
% Find lower and upper cell boundaries:
%
Delta = m(2)/m(1);
mLower = m/sqrt(Delta);
mUpper = m*sqrt(Delta);
m2 = min(max(mUpper),m2);
%
% Find affected size ranges:
%
ix = find(mUpper>=m1 & mLower<=m2);
%
% Find fractions to take from boundary size groups:
%
if mLower>m1
    fLower = 1;
else
    fLower = 1-(log(m1)-log(mLower(ix(1)))) / (log(mUpper(ix(1)))-log(mLower(ix(1))));
end

if mUpper<m2
    fUpper = 1;
else
    fUpper = (log(m2)-log(mLower(ix(end)))) / (log(mUpper(ix(end)))-log(mLower(ix(end))));
end
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
Brange = sum(B.*f); 