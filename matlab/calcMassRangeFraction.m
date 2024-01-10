%
% Calculates the fraction of each size group that is within a given range
% of masses. Only works with the size groups in the biomass range.
%
% In:
%  p - The parameter structure
%  m1, m2 - The biomass range
%
% Out:
%  The (geometric) fraction that is within the range m1 and m2 for each
%  size group.
%
function f = calcMassRangeFraction(p, m1, m2)

m1 = max(1e-20,m1); % Avoid zero lower range

f = zeros(1,p.n-p.idxB+1);

for i = p.idxB:length(p.m)
    if m2 < p.mLower(i) || m1 > p.mUpper(i)
        f(i-p.idxB+1) = 0;
    else
        mMin = max(m1, p.mLower(i));
        mMax = min(m2, p.mUpper(i));
        f(i-p.idxB+1) = (log(mMax)-log(mMin)) / (log(p.mUpper(i))-log(p.mLower(i)));
    end
    %fprintf("%f %f %f %f, %f\n",[log(mMin),log(mMax),log(p.mLower(i)),log(p.mUpper(i)),f(i-p.idxB+1)])
end
