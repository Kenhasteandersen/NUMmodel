%
% Make the numbers for table ?
%
addpath('..')

nBins = 25;
uM_to_ug = 30.97; % conversion from Molar to ugP/l.

d = [0.01 0.1]; % Mixing rates
P0 = [0.006 0.06 0.6] * uM_to_ug;

p = parametersGeneralistsOnly( nBins, 10.);
p = parametersChemostat( p );
p.mortHTLm = 0*p.mortHTLm; % No HTL mortality

p.bLosses = true; % Allow mixing losses

pp = p;

fprintf('\n');
fprintf('Total uptake of carbon (light, DOC and feeding).\n');
fprintf('\n');
fprintf('  Pico  |  Nano   |  Micro   \n');
fprintf('-----------------------------\n');
%
% Case one: no phagotrophy
%
p = pp;
p.AF = 0*p.AF; % Setting the affinity for feeding to zero
sweep(p,d,P0)
fprintf('-----------------------------------\n');
%
% Case two: pico phototrophy
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
            p.tEnd = 500;
            
            sim = simulateChemostat(p);
            plotChemostat(sim)
            drawnow
            
            m = p.m(3:end);
            Bpnm = calcPicoNanoMicro(mean(sim.B(floor(end/2):end,:),1), m);
            %
            % Average the rates over the last half of the simulation:
            %
            JLpnmIntegral = [0 0 0];
            JDOCpnmIntegral = [0 0 0];
            JFpnmIntegral = [0 0 0];
            for k = floor(length(sim.t)/2) : length(sim.t)
                u = [sim.N(k), sim.DOC(k), sim.B(k,:)];
                rates = calcDerivatives(p,u,sim.L);
                
                jL = rates.JLreal(3:end) ./m;
                jDOC = rates.JDOC(3:end) ./m;
                jF = rates.JFreal(3:end) ./m;
                %jC = rates.Jtot(3:end) ./ m;
                jLpnm = calcPicoNanoMicroRate(m, jL);
                jDOCpnm = calcPicoNanoMicroRate(m, jDOC);
                jFpnm = calcPicoNanoMicroRate(m, jF);
                Bpnm = calcPicoNanoMicro(mean(sim.B(k,:),1), m);
                JLpnmIntegral = JLpnmIntegral + jLpnm.*Bpnm*(sim.t(k)-sim.t(k-1));
                JDOCpnmIntegral = JDOCpnmIntegral + jLpnm.*Bpnm*(sim.t(k)-sim.t(k-1));
                JFpnmIntegral = JFpnmIntegral + jLpnm.*Bpnm*(sim.t(k)-sim.t(k-1));
            end
            JLpnm = JLpnmIntegral / (sim.t(end) - sim.t(floor(length(sim.t)/2)));
            JDOCpnm = JDOCpnmIntegral / (sim.t(end) - sim.t(floor(length(sim.t)/2)));
            JFpnm = JFpnmIntegral / (sim.t(end) - sim.t(floor(length(sim.t)/2)));
            %Jpnm(Jpnm<0) = 0;
            %percentages = Jpnm/sum(Jpnm)*100;
            JLpercent = JLpnm ./ sum(
            
            fprintf(' %3.1f%%  |  %3.1f%%  |  %3.1f%%   |\n', percentages);
        end
    end
end