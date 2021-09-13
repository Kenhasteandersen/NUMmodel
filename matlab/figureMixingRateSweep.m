% addpath ~/Documents/Source/NUMmodel/matlab

%
% Sweep over mixing rates and calculate gross PP
%


% The first parameters is the number of size bins. The last parameter is
% the upper size:
p = parametersGeneralistsOnly( 25, 10.);
p = parametersChemostat( p );

%p.mortHTLm = p.mortHTLm / 4; % Lowering the background mortality (not needed, actually).

M = [20 30 40]; % The thickness of the mixing layer
d = logspace(-3,log10(0.5),10);

clf
%tiledlayout(length(d),2)

% Without losses to the deep
prodGross = sweep(p,d,20) % With phagotrophy
semilogx(d, prodGross)

% With losses:
p.bLosses = true;
prodGross = sweep(p,d,20) % With phagotrophy
hold on
semilogx(d, prodGross)

legend({'Without losses','With losses of B to the deep'})
xlabel('d (day^{-1}')
ylabel('Gross PP (gC/m^2/yr)')
%p.AF = 0*p.AF; % Setting the affinity for feeding to zero
%prodGross0 = sweep(p,d,20) % without phagotrophy

%
% Sweep over mixing coeffs:
%
function prodGross = sweep(pp,d,M)


p = pp;

p.u0(3:end) = 0.0001; % for faster convergence
for i = 1:length(d)
    p.d = d(i);
    
    sim = simulateChemostat( p, 60 );
    
    jL = sim.rates.JLreal./p.m;
    prodGross(i) = M*365*1e-3 * sum(jL(3:end).*sim.B(end,:))/p.pGeneralists.epsilonL;
    
    %plotChemostat(sim)
    %drawnow
end

end