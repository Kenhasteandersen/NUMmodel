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
function [m, Bc] = calcCommunitySpectrum(p,B)

[m, idx] = sort(p.m(p.idxB:end));
Bcumm = cumsum(B(idx));
Bc = gradient(Bcumm);


% m = logspace(min(log10(p.mLower(p.idxB:end))), max(log10(p.mLower(p.idxB:end)+p.mDelta(p.idxB:end))), 100);
% Bspline = 0*m;
% B = B./log(p.mUpper(p.idxB:end)./p.mLower(p.idxB:end));
% for i = 1:p.nGroups
%     ix = p.ixStart(i):p.ixEnd(i);
%     ix = ix( B(ix-p.idxB+1)>0 );
%     
%     ixSpline = ( m>p.m(p.ixStart(i)) & m<(p.m(p.ixEnd(i))) );
% 
%     logBsplineGroup = spline(log10(p.m(ix)), log10(B(ix-p.idxB+1)), log10(m(ixSpline)));
%     Bspline(ixSpline) = Bspline(ixSpline) + 10.^logBsplineGroup;
%     loglog(m(ixSpline), 10.^logBsplineGroup)
%     hold on
% end
% 

