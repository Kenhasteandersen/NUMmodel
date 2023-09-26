function [tiles,sim] = plotGlobalBiomass(sim, sProjection)
arguments
    sim struct;
    sProjection string = 'fast';
end
%
% Calculate pico-nano-micro in each water column and average over the year (assuming monthly data):
%
if ~isfield(sim,'Bpnm')
    sim.Bpnm = zeros(length(sim.x), length(sim.y), 3);
    for iTime = find(sim.t > sim.t(end)-365) % Average over last year
        for i = 1:length(sim.x)
            for j = 1:length(sim.y)
                tmp = [0 0 0];
                for k = 1:length(sim.z)
                    tmp2 = calcPicoNanoMicro(squeeze(sim.B(iTime,i,j,k,:)), sim.p.m(sim.p.idxB:end));
                    tmp2(isnan(tmp2))=0;
                    tmp = tmp + tmp2 * sim.dznom(k)*0.001; % gC/m2
                end
                sim.Bpnm(i,j,1:3) = sim.Bpnm(i,j,1:3) + reshape(tmp,1,1,3);
            end
        end
    end
    sim.Bpnm = sim.Bpnm/iTime;
end

%%
clf
tiles = tiledlayout(3,1,'TileSpacing','compact','padding','compact');

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(sim.Bpnm(:,:,1))),[-1,1],...
    sTitle='Pico',sProjection=sProjection);
cbar.Visible='off';
caxis([-1,1])

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(sim.Bpnm(:,:,2))),[-1,1],...
    sTitle='Nano', sProjection=sProjection);
cbar.Visible='off';
caxis([-1,1])

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(sim.Bpnm(:,:,3))),[-1,1],...
    sTitle='Micro', sProjection=sProjection);
caxis([-1,1])

cbar.Label.String = 'log10(gC m^{-2})';
cbar.Location = 'SouthOutside';
cbar.FontSize = 10;



