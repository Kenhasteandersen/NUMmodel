% addpath ~/Documents/Source/NUMmodel/matlab

%
% Sweep over mixing rates and calculate gross PP
%

% The first parameters is the number of size bins. The last parameter is
% the upper size:
p = parametersGeneralistsOnly( 25, 10.);
p = parametersChemostat( p );

p.mortHTLm = p.mortHTLm/4; % Lowering the background mortality (not needed, actually).

p.u0(1) = 2;

p.tEnd = 500;

d = logspace(-3,log10(0.5),25);

clf
%tiledlayout(length(d),2)

% Without losses to the deep
prodGross1 = sweep(p,d,20) % With phagotrophy
% With losses:
p.bLosses = true;
prodGross2 = sweep(p,d,20) % With phagotrophy
%%
clf
semilogx(d, prodGross1)
hold on
semilogx(d, prodGross2)

patch([d d(end:-1:1)], [prodGross1 prodGross2(end:-1:1)],[0.8 0.8 1]);
semilogx(d, prodGross1,'b-','linewidth',2)
semilogx(d, prodGross2,'r-','linewidth',2)
ylim([0 1100])
set(gca,'xtick',10.^(-3:1))

legend({'Without losses','With losses of B to the deep'},...
    'Location','NorthWest')
legend('boxoff')
xlabel('Mixing rate (day^{-1})')
ylabel('Gross PP (gC/m^2/yr)')
xlim([min(d) max(d)])
set(gca,'fontsize',10)

% Make fig single column:
set(gcf,'paperunits','centimeters')
set(gcf,'units','centimeters')
pos = get(gcf,'paperposition');
pos(3) = 8;
set(gcf,'paperposition',pos);
ppos = get(gcf,'position');
  ppos([3 4]) = pos([3 4]);
  set(gcf,'position',ppos)
%p.AF = 0*p.AF; % Setting the affinity for feeding to zero
%prodGross0 = sweep(p,d,20) % without phagotrophy

%
% Sweep over mixing coeffs:
%
function prodGross = sweep(pp,d,M)

L = 30; % The light level to use

p = pp;

p.u0(3:end) = 0.0001; % for faster convergence
for i = 1:length(d)
    p.d = d(i);
    
    sim = simulateChemostat( p, L );
    %
    % Calc PP over the last half of the simulation:
    %
    %jL = sim.rates.JLreal./p.m;
    %prodGross(i) = M*365*1e-3 * sum(jL(3:end).*sim.B(end,:))/p.pGeneralists.epsilonL; 
    prodGross(i) = 0;
    dt = gradient(sim.t);
    ixStart = find(sim.t>(sim.t(end)/2),1);
    for j = ixStart:length(sim.t)
        u = [sim.N(j) sim.DOC(j) sim.B(j,:)];
        rates = calcDerivatives(p,u,L);
        jL = rates.JLreal./p.m;
        prodGross(i) = prodGross(i) + ...
            M*365*1e-3 * sum(jL(3:end).*sim.B(j,:))/p.pGeneralists.epsilonL * dt(j);
    end
    %plotChemostat(sim)
    %drawnow
end
prodGross = prodGross / (sim.t(end)-sim.t(ixStart));

end