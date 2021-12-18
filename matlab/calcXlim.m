%
% Returns the range of mass values:
%
function limits = calcXlim(p)

mMin = 100;
mMax = 0;

for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    mMin = min([mMin, m]);
    mMax = max([mMax, m]);
end

limits = [mMin, mMax];