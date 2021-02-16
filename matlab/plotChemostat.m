%
% Makes a series of plot of a chemostat run.
%
function plotChemostat(sim)

clf
tiledlayout(4,1,'tilespacing','compact','padding','compact')
%
% Size spectrum:
%
nexttile
panelSpectrum(sim)
%set(gca,'xticklabel','')
xlabel('')
%
% Rates:
%
nexttile;
panelGains(sim.p, sim.rates)

nexttile
panelLosses(sim.p, sim.rates)
%
% Time
%
nexttile
panelTime(sim)
