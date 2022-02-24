%
% Calculate the community size spectrum from all groups.
%
% In:
%  m - mass of all biomass groups
%  B - biomasses
%
% Out:
%  mc - masses of the community spectrum (sorted version of m)
%  Bc - Sheldon community spectrum
%
function [m, Bc] = calcCommunitySpectrum(m,B)

[m, idx] = sort(m);
Bcumm = cumsum(B(idx));
Bc = gradient(Bcumm);


% m = logspace(min(log10(p.mLower(p.idxB:end))), max(log10(p.mLower(p.idxB:end)+p.mDelta(p.idxB:end))), 100);
% Bspline = 0*m;
% B = B./log(p.mUpper(p.idxB:end)./p.mLower(p.idxB:end));
% for i = 1:p.nGroups
%     ix = p.ixStart(i):p.ixEnd(i);
%     ixSpline = ( m>p.mLower(p.ixStart(i)) & m<(p.mLower(p.ixEnd(i))+p.mDelta(p.ixEnd(i))));
% 
%     BsplineGroup = spline(p.m(ix), B(ix-p.idxB+1), m(ixSpline))
%     Bspline(ixSpline) = Bspline(ixSpline) + BsplineGroup;
%     loglog(m(ixSpline), BsplineGroup)
%     hold on
% end

