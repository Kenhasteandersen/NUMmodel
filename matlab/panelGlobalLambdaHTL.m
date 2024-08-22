%
% Plot the annual average trophic level of higher level production
% in the upper layer of a global simulation
%
function lambdaHTL = panelGlobalLambdaHTL(sim)

lambdaHTL = zeros(length(sim.t),length(sim.x),length(sim.y));
PHTL = zeros(length(sim.t),length(sim.x),length(sim.y));
nx = length(sim.x);
ny = length(sim.y);
iz = 1;

N = sim.N;
DOC = sim.DOC;
Si = sim.Si;
BB = sim.B;

parfor iTime = find(sim.t>sim.t(end)-365)
    for ix = 1:nx
        for iy = 1:ny
            if isfield(sim.p,'idxSi')
                u = [squeeze(N(iTime,ix,iy,iz)), ...
                    squeeze(DOC(iTime,ix,iy,iz)), ...
                    squeeze(Si(iTime,ix,iy,iz)), ...
                    squeeze(BB(iTime,ix,iy,iz,:))'];
            else
                u = [squeeze(N(iTime,ix,iy,iz)), ...
                    squeeze(DOC(iTime,ix,iy,iz)), ...
                    squeeze(BB(iTime,ix,iy,iz,:))'];
            end
            rates = getRates(sim.p, u, squeeze(sim.L(iTime,ix,iy,iz)), squeeze(sim.T(iTime,ix,iy,iz)));

            PHTL(iTime,ix,iy) = sum(squeeze(BB(iTime,ix,iy,iz,:)).*rates.mortHTL);
            [~, lHTL] = calcTrophicLevel(sim.p, squeeze(BB(iTime,ix,iy,iz,:)), rates);
            lambdaHTL(iTime,ix,iy) = lHTL;
        end
    end
end

% Average over time weighing with biomass:
ix = isnan(lambdaHTL);
lambdaHTL(ix) = 0;
lambdaHTL = squeeze( sum(lambdaHTL.*PHTL,1) ./ sum(PHTL,1) );

panelGlobal(sim.x, sim.y, lambdaHTL, sProjection="eckert4",sUnits='')
title('Trophic level of HTL production')


