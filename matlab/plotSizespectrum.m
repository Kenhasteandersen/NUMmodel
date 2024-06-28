%
% Plot a size spectrum at a given time (day).
% If the simulation is a watercolumn then indicate also the depth layer.
% If the simulation is global then indicate depth layer, and latitude and longitude.
% 
function s = plotSizespectrum(sim, time, iDepth, lat, lon)
arguments
    sim struct;
    time = sim.p.tEnd;
    iDepth {mustBeInteger} = 1;
    lat double = [];
    lon double = [];
end

m = [sim.p.mLower(3:end), sim.p.mLower(end)+sim.p.mDelta(end)];
[~, iTime] = min(abs(sim.t-time));

switch sim.p.nameModel
    
    case 'chemostat'
        % Extract from a chemostat:
        s.B = sim.B;
        if isfield(sim,'Si')
            u = [sim.N(iTime), sim.DOC(iTime), sim.Si(iTime), squeeze(sim.B(iTime,:))];
        else
            u = [sim.N(iTime), sim.DOC(iTime), squeeze(sim.B(iTime,:))];
        end
        s.L = mean(sim.L);
        s.T = sim.T;
        
    case 'watercolumn'
        % Extract from a single water column:
        z = sim.z;
        s.B = squeeze(sim.B(:,iDepth,:));
        if isfield(sim,'Si')
            u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), sim.Si(iTime,iDepth), ...
                squeeze(s.B(iTime,:))];
        else
            u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), squeeze(s.B(iTime,:))];
        end
        s.L = sim.L(iTime,iDepth);
        s.T = sim.T(iTime,iDepth);
        
    case 'global'
        % Extract from global run:
        idx = calcGlobalWatercolumn(lat,lon,sim);
        z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];
        s.B = squeeze(sim.B(:, idx.x, idx.y, iDepth, :));
        if isfield(sim,'Si')
            u = [sim.N(iTime,idx.x, idx.y,iDepth), ...
                sim.DOC(iTime,idx.x, idx.y,iDepth), ...
                sim.Si(iTime,idx.x, idx.y,iDepth), ...
                squeeze(sim.B(iTime,idx.x, idx.y, iDepth, :))'];
        else
            u = [sim.N(iTime,idx.x, idx.y,iDepth), ...
                sim.DOC(iTime,idx.x, idx.y,iDepth), ...
                squeeze(sim.B(iTime,idx.x, idx.y, iDepth, :))'];
        end
        s.L = sim.L(iTime,idx.x, idx.y, iDepth);
        s.T = sim.T(iTime,idx.x, idx.y, iDepth);
end

s.p = sim.p;
s.t = sim.t;
if ~isfield('sim','rates')
    sim.rates = getRates(sim.p, u, s.L, s.T);
end
s.rates = sim.rates;
%
% Setup tiles:
%
clf
tiledlayout(4,1,'tilespacing','compact','padding','compact')
%
% Spectrum
%
nexttile
panelSpectrum(s,iTime,bPlotStrategies=false)
xlabel('')
set(gca,'XTickLabel','');
%
% Gains:
%
nexttile

%panelGains(sim.p,rates)
panelGains(sim.p, sim.rates);
set(gca,'XTickLabel','');
%
% Losses:
%
nexttile
panelLosses(sim.p, sim.rates);
set(gca,'XTickLabel','');
xlabel('')

nexttile
panelTrophicLevel(sim,sim.rates,lat,lon);


if strcmp(sim.p.nameModel, 'watercolumn')

sgtitle(['Day = ', num2str(time), ', lat = ', num2str(sim.lat), char(176), ', lon = ', num2str(sim.lon), char(176), ', depth of maximum biomass: ', num2str(z(iDepth)), ' m']) 

end