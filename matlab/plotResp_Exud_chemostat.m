 p = setupGeneralistsOnly;
p = setupGeneralistsDiatoms;
% p = setupGenDiatCope(2,2,1);
%  p = setupDiatomsOnly;
p = parametersChemostat(p);
p.tEnd = 2000;
% p.d = 0.1;
%
% Set to "normal" HTL mortality if there are no copepods:
% %
% if isempty(mAdult)
%     setHTL(0.1, 1/500^1.5, false, false);
% else 
    setHTL(0.1, 1, true, true);
% end
%%
% Simulate
%
%-----------    eutrophic   -------------
p.d=0.1;
tic
simEut = simulateChemostat(p, 40);
toc
%
%% Plot
%
group=1;
plotSimulation(simEut);
figure(2)
% nexttile
panelRespirationExudation(simEut.p,simEut.rates, group)
% nexttile
% panelExudation(simEut.p,simEut.rates, group)
%%
%-----------    oligotrophic   -------------
p.d=0.001;
tic
simOl = simulateChemostat(p, 100);
toc
%%
group=1;

%
% Plot
%
plotSimulation(simOl);
figure(3)
% nexttile
panelRespirationExudation(simOl.p,simOl.rates, group)

% time = datestr(clock,'YYYY_mm_dd_HH_MM_SS');
% % 
% saveas(gcf,['panelResp_',time,'.png'])