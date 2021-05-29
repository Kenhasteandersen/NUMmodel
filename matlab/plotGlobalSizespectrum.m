%
% Plot a size spectrum at a given lat, lon, depth (in meters), and time
% (index).
%
function s = plotGlobalSizespectrum(sim, lat, lon, depth, iTime)
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

if isfield(sim,'Si')
u = [sim.N(idx.x, idx.y,iDepth,iTime), ...
    sim.DOC(idx.x, idx.y,iDepth,iTime), ...
    sim.Si(idx.x, idx.y,iDepth,iTime), ...
    squeeze(sim.B(idx.x, idx.y, iDepth, :, iTime))'];
else
    u = [sim.N(idx.x, idx.y,iDepth,iTime), ...
        sim.DOC(idx.x, idx.y,iDepth,iTime), ...
        squeeze(sim.B(idx.x, idx.y, iDepth, :, iTime))'];
end

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

s.L = sim.L(idx.x, idx.y, iDepth,iTime)
rates = getRates(sim.p, u, s.L);
panelGains(sim.p,rates)
%
% Losses:
%
nexttile
panelLosses(sim.p, rates);

s.rates = rates;