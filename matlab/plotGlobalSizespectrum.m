%
% Plot a size spectrum at a given lat, lon, and depth
%
function plotGlobalSizespectrum(sim, lat, lon, depth, iTime)
arguments
    sim struct;
    lat, lon, depth double;
    iTime {mustBeInteger} = length(sim.t);
end

idx = calcGlobalWatercolumn(lat,lon,sim);
m = [sim.p.mLower(3:end), sim.p.mLower(end)+sim.p.mDelta(end)];
z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];
[~, iDepth] = min(abs(z-depth));

s.p = sim.p;
s.B = squeeze(sim.B(idx.x, idx.y, iDepth, :, iTime))';
s.t = sim.t;


clf
tiledlayout(3,1,'tilespacing','compact','padding','compact')
%
% Spectrum
%
nexttile
panelSpectrum(s,1)

%
% Gains:
%
nexttile
u = [sim.N(idx.x, idx.y,iDepth,iTime), ...
    sim.DOC(idx.x, idx.y,iDepth,iTime), ...
    squeeze(sim.B(idx.x, idx.y, iDepth, :, iTime))'];
rates = getRates(u, sim.L(idx.x, idx.y, iDepth,iTime));
panelGains(sim.p,rates)
%
% Losses:
%
nexttile
panelLosses(sim.p, rates);