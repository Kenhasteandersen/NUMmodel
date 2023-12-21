%
% Calculate the community size spectrum from all groups using interplotation.
%
% In:
%  B - vector of biomasses (time, biomass)
%  sim - simulation structure
%  iTime - (optional) time index to plot (if not set an average is used)
%
% Out:
%  mc - central masses of bins
%  BSheldon - Sheldon size spectrum as defined in Andersen and Visser (2023), Box V.
%  Bspectrum - Normalized biomass spectrum
%
function [mc, BSheldon, Bspectrum] = calcCommunitySpectrum(B, sim, iTime)

arguments
    B;
    sim struct;
    iTime = NaN;
end

B(B<=0) = 1e-100; % just to avoid imaginary numbers during log transformation

p = sim.p;

nPoints = 1000;
mc = logspace(log10(min(sim.p.m(p.idxB:end))), log10(max(sim.p.m)), nPoints);
BSheldon = zeros(1, nPoints);

for iGroup = 1:p.nGroups

    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    Delta = p.mUpper(ix)./p.mLower(ix);
    ixB = ix-p.idxB+1;

    if isnan(iTime)
        ixAve = find( sim.t > sim.t(end)/2 );
        % Interpolation
        log_k1 = mean( log(B( ixAve, ixB)./log(Delta)),1);
        if length(ix)==1
            vq1 = exp(interp1( log([p.mLower(ix) p.mUpper(ix)]), [log_k1, log_k1], log(mc),'linear') );
        else
            vq1 = exp(interp1(log(m), log_k1, log(mc), 'linear'));
        end

        vq1(isnan(vq1)) = 0; % get rid of the NAs
        BSheldon = BSheldon + vq1;
    else

        % Interpolation
        log_k1 = mean( log(B( iTime, ixB)./log(Delta)),1);
        if length(m)==1
            vq1 = exp(interp1( log([p.mLower(ix) p.mUpper(ix)]), [log_k1, log_k1], log(mc), 'linear'));
        else
            vq1 = exp(interp1(log(m), log_k1, log(mc), 'linear'));
        end

        vq1(isnan(vq1)) = 0; % get rid of the NAs
        BSheldon = BSheldon + vq1;
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

Delta = mc(2)./mc(1);
dff = diff(log(mc));
mUpper = exp( log(mc(1:end-1)) + 0.5*dff);
mUpper(end+1) = exp( log(mc(end)) + 0.5*dff(end));
mLower(2:length(mUpper)) = mUpper(1:end-1);
mLower(1) = exp( log(mc(1)) - 0.5*dff(1));
Bspectrum = BSheldon ./ (mUpper-mLower) * Delta;

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