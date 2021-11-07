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
fprintf('         Pico         |         Nano          |         Micro          \n');
fprintf(' Light    DOC   Food  |  Light   DOC    Food  | Light   DOC    Food    \n');
fprintf('-----------------------------------------------------------------------\n');
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
            p.tEnd = 1000;
            
            sim = simulateChemostat(p);
            plotChemostat(sim)
            drawnow
            
            m = p.m(3:end);
            %
            % Average the rates over the last half of the simulation:
            %
            JLpnmIntegral = [0 0 0];
            JDOCpnmIntegral = [0 0 0];
            JFpnmIntegral = [0 0 0];
            for k = find(sim.t>0.75*sim.t(end),1):length(sim.t)
                u = [sim.N(k), sim.DOC(k), sim.B(k,:)];
                rates = calcDerivatives(p,u,sim.L);
                % Get the rates:
                jL = rates.JLreal(3:end) ./m;
                jDOC = rates.JDOC(3:end) ./m;
                jF = rates.JFreal(3:end) ./m;
                jC = jL + jDOC + jF;
                % Sum them in pico-nano-micro:
                jLpnm = calcPicoNanoMicroRate(m, jL);
                jDOCpnm = calcPicoNanoMicroRate(m, jDOC);
                jFpnm = calcPicoNanoMicroRate(m, jF);
                % Calc the percentages:
                Bpnm = calcPicoNanoMicro(sim.B(k,:), m);
                JLpnmIntegral = JLpnmIntegral + jLpnm.*Bpnm;
                JDOCpnmIntegral = JDOCpnmIntegral + jDOCpnm.*Bpnm;
                JFpnmIntegral = JFpnmIntegral + jFpnm.*Bpnm;
                %jCpnmIntegration = JCpnmIntegration + jC.*Bpnm*deltaT;
            end
            %deltaT = 1;%sim.t(end) - sim.t(floor(length(sim.t)/2));
            JC = sum(JLpnmIntegral + JDOCpnmIntegral + JFpnmIntegral);
            JLpnm = 100*JLpnmIntegral/JC;
            JDOCpnm = 100*JDOCpnmIntegral/JC;
            JFpnm = 100*JFpnmIntegral/JC;
                        
            fprintf(' % 3.1f%% % 3.1f%% % 3.1f%%  |  % 3.1f%% % 3.1f%% % 3.1f%% |  % 3.2f%% % 3.2f%% % 3.2f%%  | %3.1f%%\n', ...
                [JLpnm(1), JDOCpnm(1), JFpnm(1), ...
                JLpnm(2), JDOCpnm(2), JFpnm(2), ...
                JLpnm(3), JDOCpnm(3), JFpnm(3), ...
                sum(JLpnm+JDOCpnm+JFpnm)]);
            
            %ix = sim.t>(0.75*sim.t(end));
            %Bpnm = calcPicoNanoMicro(trapz(sim.t(ix),sim.B(ix,:),1), p.m(3:end)) / (0.25*sim.t(end));
            %fprintf(' %3.1f  |  %3.1f  |  %3.1f   |\n', Bpnm);
        end
    end
end