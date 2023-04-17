function [mc, Bc] = calcCommunitySpectrumWaterCol(B, sim,time)
% θ
arguments
    B;
    sim struct;
    time = NaN;
end

p = sim.p;

nPoints = 1000;
mc = logspace(log10(sim.p.m(p.idxB)), log10(sim.p.m(end)), nPoints);
Bc = zeros(1, nPoints);
ixAve = find( sim.t > sim.t(end)/2 );

for iGroup = 1:p.nGroups

    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    Delta = p.mUpper(ix)./p.mLower(ix);
    ixB = ix-p.idxB+1;


    if isnan(time)
        % Interpolation
        log_k = mean( log(B(ixAve, ixB)./log(Delta)),1);
        vq1 = exp(interp1(log(m), log_k, log(mc), 'linear'));

        vq1(isnan(vq1)) = 0; % get rid of the NAs
        Bc = Bc + vq1;

    else

        log_k = mean( log(B(ixB)./log(Delta)),1);
        vq1 = exp(interp1(log(m), log_k, log(mc), 'linear'));

        vq1(isnan(vq1)) = 0; % get rid of the NAs
        Bc = Bc + vq1;


    end
end
end


% function Bc = calcCommunitySpectrum(sim, mc)
%
%     p = sim.p;
%     nPoints = length(mc);
%     Bc = zeros(1, nPoints);
%
%     for iGroup = 1:p.nGroups
%
%         ix = p.ixStart(iGroup):p.ixEnd(iGroup);
%         m = p.m(ix);
%         Delta = p.mUpper(ix)./p.mLower(ix);
%         ixB = ix-p.idxB+1;
%         ixAve = find( sim.t > sim.t(end)/2 );
%
%         % Interpolation
%         log_k = mean( log(sim.B(ixAve, ixB)./log(Delta)),1);
%         vq1 = exp(interp1(log(m), log_k, log(mc), 'linear'));
%
%         vq1(isnan(vq1)) = 0; % get rid of the NAs
%         Bc = Bc + vq1;
%     end
%

