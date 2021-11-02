%
% Make the numbers for table one
%
nBins = 25;
uM_to_ug = 30.97; % conversion from Molar to ugP/l.

d = [0.01 0.1]; % Mixing rates
P0 = [0.006 0.06 0.6] * uM_to_ug;

p = parametersGeneralistsOnly( nBins, 10.);
p = parametersChemostat( p );
p.mortHTLm = 0*p.mortHTLm; % No HTL mortality

p.bLosses = true; % Allow mixing losses

pp = p;

fprintf('  Pico  |  Nano   |  Micro |  ugC/L\n');
fprintf('-----------------------------------\n');
%
% Case one: no phagotrophy
%
p = pp;
p.AF = 0*p.AF; % Setting the affinity for feeding to zero
sweep(p,d,P0)
fprintf('-----------------------------------\n');
%
% Case two: only pico phototrophy
%
p = pp;
p.pGeneralists.ALm( p.m(3:end)>1e-6 ) = 0; % No phototrophy for larger cells
sweep(p,d,P0)
fprintf('-----------------------------------\n');
%
% Case three: full model
%
p = pp;
sweep(p,d,P0)
fprintf('-----------------------------------\n');

function sweep(p,d,P0)
    for i = 1:length(d)
        for j = 1:length(P0)
            p.d = d(i);
            p.u0(1) = P0(j);
            p.tEnd = 1000;
            
            sim = simulateChemostat(p);
            plotChemostat(sim)
            drawnow
            
            ix = sim.t>(0.75*sim.t(end));
            Bpnm = calcPicoNanoMicro(mean(sim.B(ix,:),1), p.m(3:end));
            fprintf(' %3.1f  |  %3.1f  |  %3.1f   |\n', Bpnm);
        end
    end
end