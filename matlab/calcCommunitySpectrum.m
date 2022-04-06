%
% Calculate the community size spectrum from all groups.
%

function [mc, Bc] = calcCommunitySpectrum(sim)

p = sim.p;
mX = p.m;
nPoints = 100000;
mc = logspace(log10(mX(3)), log10(mX(end)), nPoints);
    
Bc = zeros(1, nPoints);
    
for iGroup = 1:p.nGroups
    
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    Delta = p.mUpper(ix)./p.mLower(ix);
    ixB = ix-p.idxB+1;
    ixAve = find( sim.t > sim.t(end)/2 );
    
    % Interpolation
    log_k = mean( log(sim.B(ixAve, ixB)./log(Delta)),1);
    vq1 = exp(interp1(log(m), log_k, log(mc), 'linear'));
    
    vq1(isnan(vq1)) = 0;
    Bc = Bc + vq1;
    end

end



% [m, idx] = sort(p.m(3:end));
% Bcumm = cumsum(B(idx));
% Bc = gradient(Bcumm);

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

