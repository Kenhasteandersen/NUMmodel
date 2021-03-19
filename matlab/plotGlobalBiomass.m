function plotGlobalBiomass(sim, sProjection)
arguments
    sim struct;
    sProjection string = 'fast';
end
%
% Calculate pico-nano-micro in each water column and average over the year (assuming monthly data):
%
Bpnm = zeros(length(sim.x), length(sim.y), 3);
for iTime = 1:12
    for i = 1:length(sim.x)
        for j = 1:length(sim.y)
            tmp = [0 0 0];
            for k = 1:length(sim.z)
                tmp2 = calcPicoNanoMicro(squeeze(sim.B(i,j,k,:,iTime)), sim.p);
                tmp2(isnan(tmp2))=0;
                tmp = tmp + tmp2;
            end
            Bpnm(i,j,1:3) = Bpnm(i,j,1:3) + reshape(tmp/12,1,1,3);
        end
    end
end
%%
clf
tiledlayout(1,3,'TileSpacing','compact','padding','compact')

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(Bpnm(:,:,1))),'Pico',sProjection);
cbar.Visible='off';
caxis([-3,2])

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(Bpnm(:,:,2))),'Nano',sProjection);
cbar.Visible='off';
caxis([-3,2])

nexttile
panelGlobal(sim.x,sim.y,(log10(Bpnm(:,:,3))),'Micro',sProjection);
caxis([-3,2])
cbar.Label.String  = 'g C m{^-2};'

