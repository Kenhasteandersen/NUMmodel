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

s.p = sim.p;
s.t = sim.t;
if ~isfield('sim','rates')
    sim.rates = getRates(sim.p, u, s.L, s.T);
end
s.rates = sim.rates;


[Bphyto, Bzoo, Bbacteria] = calcPhytoZoo(s.p, u, s.L, s.T);
%%
figure(33)
clf
plot(Bphyto,Bzoo,'o')

%%
t=linspace(1,10,50);
newPalette={[21, 16, 240]/256,[0, 163, 136]/256,[110, 85, 71]/256,...
    [219, 6, 6]/256,[219, 6, 6]/256,[219, 6, 6]/256,[219, 6, 6]/256,[116, 71, 145]/256};
figure(100)
clf
tiledlayout(1,2)
nexttile
for i=1:8
    hold on
 plot(t,t+i,'LineWidth',2,'Color',newPalette{i})
end

nexttile
for i=1:8
    hold on
 plot(t,t+i,'LineWidth',2,'Color',p.colGroup{i})
end
%%
% nesting tiles

T=tiledlayout(3,1); %Outer layout
nexttile(T,[1,1]); axis off %Next outer tile
t=tiledlayout(T,1,2); %first inner layout
for i=1:2
  nexttile(t);
  plot(rand(5,m1(i)));
  legend('Location','southoutside')
end
nexttile(T,[2,1]); axis off   %Next outer tile
t=tiledlayout(T,2,2); %second inner layout
t = tiledlayout(2,2);
t.Subtitle.String = 'An Insightful Subtitle';
t.Subtitle.FontAngle = 'italic';