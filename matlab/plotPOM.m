% 
% Plot basic POM plots - Author: Cécile Decker
%
% Warning: the "setup" function as well and a certain simulation and the 
% calcPOM function need to be called before this plot is made. If not 
% matlab may crash or the results be incorrect.
%
% In:
%  sim - a certin simulation with a set up containing POM
%  time - time (in days)
%

function plotPOM(sim, option)

arguments
    sim;
    option.time double = 300;
end

POM = sim.POM;

%figure(3)
%plot(log10(POM.r),POM.distrib);
%figure(4)
%plot(log10(POM.r), POM.w);

figure(1)
t = sim.t;
% Make a layer at z = 0 with the same value as in the first grid point:
z = POM.Z;
N = POM.JN;
DOC = POM.JDOC;
DIC = POM.JDIC;
B = squeeze(sum(POM.B_i,3));
%ylimit = [-options.depthMax, 0];
xlimit = [sim.t(1) sim.t(end)];
lat = 60;
lon = -10;
% Nitrogen:
nexttile
contourf(t,-z,log10(N).','LineStyle','none') %linspace(-2,3,options.nLevels),
title('Nitrogen flux')
ylabel('Depth (m)')
axis tight
colorbar %('ticks',-2:3)
%ylim(ylimit)
xlim(xlimit)
xlabel('Time (days)')
set(gca,'XTickLabel','')
% DOC
nexttile
contourf(t,-z,log10(DOC).','LineStyle','none') %,options.nLevels
title('DOC flux')
ylabel('Depth (m)')
shading interp
axis tight
colorbar
%ylim(ylimit)
xlim(xlimit)
xlabel('Time (days)')
set(gca,'XTickLabel','')
% DOC
nexttile
contourf(t,-z,log10(DIC).','LineStyle','none') %,options.nLevels
title('DIC flux')
ylabel('Depth (m)')
shading interp
axis tight
colorbar
%ylim(ylimit)
xlim(xlimit)
xlabel('Time (days)')
set(gca,'XTickLabel','')
% POM
nexttile
%i = sim.p.ixPOM;
B(B < 0.01) = 0.01; % Set low biomasses to the lower limit to avoid white space in plot
contourf(t,-z,log10(B).','LineStyle','none') %,[linspace(-2,3,options.nLevels)]
title('POM Biomass')
ylabel('Depth (m)')
axis tight
colorbar
%colorbar('ticks',-12:-2) %'limits',[-12 -2]
%ylim(ylimit)
xlim(xlimit)
xlabel('Time (days)')

if strcmp(sim.p.nameModel, 'watercolumn')
    sgtitle(['Nitrogen, DOC & DIC fluxes ({\mu}gC/l/day) from POM/ Total biomass of POM ({\mu}gC/l) at: lat = ', num2str(lat), char(176), ', lon = ', num2str(lon), char(176)])
end

%% day 600 biomass as a function of size

figure(2)
time = option.time;
nexttile
[~, iTime] = min(abs(sim.t-time));
% Extract data from water column simulation:
z = POM.Z;
B = squeeze(POM.B_i(iTime,:, :));
% Biomass spectrum:
contourf( POM.r, -z, B, 'linestyle','none')
set(gca,'xscale','log','colorscale','log')
xlabel("size \mu m")
%caxis([0.01 100])
ylabel('Depth (m)')
%ylim(ylimit)
cbar = colorbar;
cbar.Label.String  = 'POM biomass (\mug C l^{-1})';

if strcmp(sim.p.nameModel,'watercolumn') || strcmp(sim.p.nameModel,'global')
    sgtitle(['Day: ', num2str(time), ', lat = ', num2str(sim.lat), char(176), ', lon = ', num2str(sim.lon), char(176)])
else
    sgtitle(['Day: ', num2str(time)])
end

%{
nexttile
% Create Sheldon biomass spectrum:
m = POM.m;
m_range(1) = m(1);
for i=1:(length(m)-1)
    m_range(i+1) = (m(i+1)-m(i))/2;
    Delta(i) = m_range(i+1)/m_range(i);
end
Delta(length(m)) = (m(length(m))+ m_range(length(m)))/m_range(length(m));
B_sheldon = B./log(Delta);
contourf( m, -z, B_sheldon, 'linestyle','none')
set(gca,'xscale','log','colorscale','log')
xlabel("mass")
%caxis([0.01 100])
ylabel('Depth (m)')
ylim(ylimit)
cbar = colorbar;
cbar.Label.String  = 'POM "Sheldon" biomass (\mug C l^{-1})';

sgtitle(['Day: ', num2str(time), ', lat = ', num2str(sim.lat), char(176), ', lon = ', num2str(sim.lon), char(176)])
%}
