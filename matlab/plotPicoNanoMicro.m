% plot pico nano micro

rGroups=calcRadiusGroups(sim.p);
r_pico=find(rGroups<1);
Bpico=sim.B(:,:,:,:,r_pico);
Bpico_total=squeeze(sum(Bpico,5));
%%
figure(1)
clf
tiledlayout(3,1)

nexttile
cbar = panelGlobal(sim.x,sim.y, squeeze(log10(mean(Bpico_total,1)/1000)), [-1, 1], sProjection='Eckert4');
cbar.Label.String = "g_C/m^2";
title('pico')
%%
nexttile
cbar = panelGlobal(sim.x,sim.y, squeeze(log10(mean(sim.Bnano,1)/1000)), [-1, 1], sProjection='Eckert4');
cbar.Label.String = "g_C/m^2";
title('nano')

nexttile
cbar = panelGlobal(sim.x,sim.y, squeeze(log10(mean(sim.Bmicro,1)/1000)), [-1, 1], sProjection='Eckert4');
cbar.Label.String = "g_C/m^2";
title('micro')