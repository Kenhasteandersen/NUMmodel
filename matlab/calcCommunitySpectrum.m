%
% Calculate the community size spectrum from all groups using interplotation.
%
function [mc, Bc] = calcCommunitySpectrum(B, sim, iTime)

arguments
    B;
    sim struct;
    iTime = NaN;
end

B(B<=0) = 1e-100; % just to avoid imaginary numbers during log transformation

p = sim.p;

nPoints = 1000;
mc = logspace(log10(min(sim.p.m(p.idxB:end))), log10(max(sim.p.m)), nPoints);
Bc = zeros(1, nPoints);

for iGroup = 1:p.nGroups

    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    Delta = p.mUpper(ix)./p.mLower(ix);
    ixB = ix-p.idxB+1;

    if isnan(iTime)
        ixAve = find( sim.t > sim.t(end)/2 );
        % Interpolation
        log_k = mean( log(B( ixAve, ixB)./log(Delta)),1);
        if length(ix)==1
            vq1 = exp(interp1( log([p.mLower(ix) p.mUpper(ix)]), [log_k, log_k], log(mc),'linear') );
        else
            vq1 = exp(interp1(log(m), log_k, log(mc), 'linear'));
        end

        vq1(isnan(vq1)) = 0; % get rid of the NAs
        Bc = Bc + vq1;
    else

        % Interpolation
        log_k = mean( log(B( iTime, ixB)./log(Delta)),1);
        if length(m)==1
            vq1 = exp(interp1( log([p.mLower(ix) p.mUpper(ix)]), [log_k, log_k], log(mc), 'linear'));
        else
            vq1 = exp(interp1(log(m), log_k, log(mc), 'linear'));
        end

        vq1(isnan(vq1)) = 0; % get rid of the NAs
        Bc = Bc + vq1;
    end

% Check interpolation (SOS! Only with marker!!!)
% figure
% 
% loglog(m, exp(log_k), '*')
% hold on 
% 
% loglog(mc, vq1)
% xlim([0, max(m)])


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
% end