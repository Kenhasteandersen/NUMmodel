%
% Figure showing selected spectra from the nutrient sweep
%
p = parametersGeneralistsOnly( 25, 10.);
p = parametersChemostat( p );
L = 30; % Light level

p.bLosses = false; % Whether to employ mixing losses on big plankton or not

p.mortHTLm = p.mortHTLm / 4; % Lowering the background mortality (not needed, actually).
p.d = 0.005; % Mixing rate

N0 = [0.1 1 10]; % The deep nutrient levels to use

AF = p.AF;
clf
for i = 1:length(N0)
   p.u0(1) = N0(i);
   
   p.AF = AF; % With phagotrophy
   sim = simulateChemostat(p,L);
   
   p.AF = 0*AF; % No phagotrophy
   simNoPhago = simulateChemostat(p,L);
   
   m = p.m(3:end);
   loglog(m, sim.B(end,:),'r-', 'linewidth',i);
   hold on
   loglog(m, simNoPhago.B(end,:),'b-', 'linewidth',i);
end

ylim([0.1 1000])
legend({'With phagotrophy','Without phagotrophy'})
xlabel('Mass ({\mu}g)')
ylabel('Biomass MUST BE NORMALIZED')
