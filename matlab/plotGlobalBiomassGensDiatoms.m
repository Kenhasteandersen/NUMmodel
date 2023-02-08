function [tiles,sim] = plotGlobalBiomassGensDiatoms(sim, sProjection)
arguments
    sim struct;
    sProjection string = 'fast';
end
%
% Calculate pico-nano-micro in each water column and average over the year (assuming monthly data):
%
if length(sim.p.ixStart)>=2        % we have at least 2 groups
    sim.Bgen=sim.B(sim.p.ixStart(1):sim.p.ixEnd(1));
    sim.Bdiat=sim.B(sim.p.ixStart(2):sim.p.ixEnd(2));
    sim.Bcope=sim.B(sim.p.ixStart(3):sim.p.ixEnd(end));
end 
if ~isfield(sim,'Bpnm')
    sim.Bpnm = zeros(length(sim.x), length(sim.y), 3);
    sim.BpnmGen=sim.Bpnm;
   if any(sim.p.typeGroups(:)==3)
     sim.BpnmDiat=sim.Bpnm;
   end
   if any(sim.p.typeGroups(:)==10)
     sim.BpnmCope=sim.Bpnm;
   end
        for iTime = find(sim.t > sim.t(end)-365) % Average over last year
            for i = 1:length(sim.x)
                for j = 1:length(sim.y)
                     tmp = [0 0 0];
          %here we need to distinguish gens, diatoms and copepods
                    for k = 1:length(sim.z)
                      sim.Bgen=sim.B(sim.p.ixStart(1):sim.p.ixEnd(1));
                      tmp2 = calcPicoNanoMicro(squeeze(sim.Bgen(i,j,k,:,iTime)), sim.p.m(sim.p.ixStart(1):sim.p.ixEnd(1)));
                      tmp2(isnan(tmp2))=0;
                      tmp = tmp + tmp2 * sim.dznom(k)*0.001; % gC/m2
                      if any(sim.p.typeGroups(:)==3)
                          sim.Bdiat=sim.B(sim.p.ixStart(2):sim.p.ixEnd(2));
                          tmp2 = calcPicoNanoMicro(squeeze(sim.Bdiat(i,j,k,:,iTime)), sim.p.m(sim.p.ixStart(2):sim.p.ixEnd(2)));
                          tmp2(isnan(tmp2))=0;
                          tmp = tmp + tmp2 * sim.dznom(k)*0.001; % gC/m2                      end
                      end
                    end
                       sim.BpnmGen(i,j,1:3) = sim.BpnmGen(i,j,1:3) + reshape(tmp,1,1,3);
                       sim.BpnmDiat(i,j,1:3) = sim.BpnmDiat(i,j,1:3) + reshape(tmp,1,1,3);
               end
            end
    sim.BpnmGen = sim.BpnmGen/iTime;
    sim.BpnmDiat = sim.BpnmDiat/iTime;
    end
end
%%
clf
tiles = tiledlayout(3,1,'TileSpacing','compact','padding','compact');

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(sim.BpnmGen(:,:,1))),[-1,1],...
    sTitle='Pico',sProjection=sProjection);
cbar.Visible='off';
caxis([-1,2])

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(sim.BpnmGen(:,:,2))),[-1,1],...
    sTitle='Nano', sProjection=sProjection);
cbar.Visible='off';
caxis([-1,2])

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(sim.BpnmGen(:,:,3))),[-1,1],...
    sTitle='Micro', sProjection=sProjection);
caxis([-1,2])

cbar.Label.String = 'log10(gC m^{-2})';
cbar.Location = 'SouthOutside';
cbar.FontSize = 10;

if any(sim.p.typeGroups(:)==3)
 figure
 tiles = tiledlayout(3,1,'TileSpacing','compact','padding','compact');

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(sim.BpnmDiat(:,:,1))),[-1,1],...
    sTitle='Pico diatoms',sProjection=sProjection);
cbar.Visible='off';
caxis([-1,2])

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(sim.BpnmDiat(:,:,2))),[-1,1],...
    sTitle='Nano diatoms', sProjection=sProjection);
cbar.Visible='off';
caxis([-1,2])

nexttile
cbar = panelGlobal(sim.x,sim.y,(log10(sim.BpnmDiat(:,:,3))),[-1,1],...
    sTitle='Micro diatoms', sProjection=sProjection);
caxis([-1,2])

cbar.Label.String = 'log10(gC m^{-2})';
cbar.Location = 'SouthOutside';
cbar.FontSize = 10;
end

end


