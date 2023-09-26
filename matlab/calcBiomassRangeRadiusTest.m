%
% Return the biomass within the size range r1:r2.
% Also output the fraction of biomass from each size
% group that falls into the range.
%
function [Brange,f] = calcBiomassRangeRadiusTest(B,r, r1, r2)
B = reshape(B,1,length(B));

%
% Find lower and upper cell boundaries:
%
Delta = r(2)/r(1);
rLower = r/sqrt(Delta);
rUpper = r*sqrt(Delta);
r2 = min(max(rUpper),r2);
%
% Find affected size ranges:
%
ix = find(rUpper>=r1 & rLower<=r2);
%
% Find fractions to take from boundary size groups:
%
if rLower>r1
    fLower = 1;
else
    if ~isempty(ix)
     fLower = 1-(log(r1)-log(rLower(ix(1)))) / (log(rUpper(ix(1)))-log(rLower(ix(1))));
    end
end

if rUpper<r2
    fUpper = 1;
else
   if ~isempty(ix)
    fUpper = (log(r2)-log(rLower(ix(end)))) / (log(rUpper(ix(end)))-log(rLower(ix(end))));
   end
end
%
% Assemble fractions:
%
f = zeros(1,length(r));
if ~isempty(ix)
 f(ix) = 1;
 f(ix(1)) = fLower;
 f(ix(end)) = fUpper;
end
%
% Assemble biomass
%
Brange = sum(B.*f); 