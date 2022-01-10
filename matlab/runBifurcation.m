%
% Creates a bifurcation plot over a field in the parameter structure.
% The parameters needs to be set up for parallel execution (set second
% arguments in call to "setupXXX" to true).
%
% In:
%  p - Parameters
%  sField - String with the field to do the bifurcation over
%  range - The value to do the bifurcation over.
%
% Example:
%  p = parametersChemostat( setupGeneric([0.1, 100], true));
%  runBifurcation(p,'d',logspace(-4,-1,10));
%
function runBifurcation(p, sField, range)

arguments
    p struct;
    sField char = 'd';
    range double = logscale(-4, -1, 10);
end

L = 30;
T = 10;
n = length(range);
%
% Set up the parameter fields:
%
for i = 1:n
    pp(i) = p;
    eval( strcat('pp(i).',sField, '= range(i);') );
end
%
% Do the simulations:
%
B = zeros(n, p.nGroups);
Blower = B;
Bupper = B;
nGroups = p.nGroups;
parfor i = 1:n
    %
    % Set up simulation:
    %
    ppp = pp(i);
    ppp.tEnd = 1000;
    %
    % Simulate:
    %
    sim = simulateChemostat(ppp, L, T);
    %
    % Calc statistics:
    %
    ixAve = find(sim.t > sim.t(end)/2);
    
    N(i) = exp( mean( log(sim.N(ixAve))));
    Nlower(i) = min(sim.N(ixAve));
    Nupper(i) = max(sim.N(ixAve));
    
    for iGroup = 1:nGroups
        ix = (ppp.ixStart(iGroup):ppp.ixEnd(iGroup)) - ppp.idxB + 1;
        B(i,iGroup) = exp(mean(log(sum(sim.B(ixAve,ix),2))));
        Blower(i,iGroup) = min(sum( sim.B(ixAve,ix),2 ));
        Bupper(i,iGroup) = max(sum( sim.B(ixAve,ix),2 ));
    end
    
    %Bpnm(i,:) = calcPicoNanoMicro(sim);
end
%%
% Plot
%

%
% Nutrients:
%
semilogy(range, N,'linewidth',2,'color','b');
hold on
patch([range range(end:-1:1)], [Nlower, Nupper(end:-1:1)], 0.75*[0,0,1], ...
    'edgecolor','none','facealpha',0.15)

legendentries(1) = semilogy(range, N,'linewidth',2,'color','b');
sLegend{1} = 'N';
%
% Biomass:
%
for iGroup = 1:nGroups
    patch([range range(end:-1:1)], [Blower(:,iGroup)', Bupper(end:-1:1,iGroup)'],...
        p.colGroup{iGroup}, 'edgecolor','none','facealpha',0.15)
   
    lwd = 2;
    if (p.typeGroups(iGroup) == 10)
        lwd = 1 + (iGroup-2)*0.5;
    end
    legendentries(1+iGroup) = ...
        semilogy(range, B(:,iGroup),'color', p.colGroup{iGroup},'linewidth',lwd);
    
    sLegend{iGroup+1} = p.nameGroup{iGroup};
end
set(gca,'xscale','log','yscale','log')
hold off
%
% Legend
%
legend(legendentries, sLegend, 'location','northwest','box','off')
xlabel(sField)
ylabel('Biomass ({\mu}g/l)')
ylim([1e-10 1e10])