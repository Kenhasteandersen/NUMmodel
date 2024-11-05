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
% Out:
%  mc, Bc - community spectra from calls to calcCommunitySpectrum
%
% Example:
%  p = parametersChemostat( setupGeneric([0.1, 100], true));
%  runBifurcation(p,'d',logspace(-4,-1,10));
%
function [mc, Bc] = runBifurcation(p, sField, range, options)

arguments
    p struct;
    sField char = 'd';
    range double = logspace(-5, 0, 10);
    options.bParallel = false;
    options.bPlotSpectra = true; % Plot also the community size spectrum
    % These options are only for the HTL bifurcation:
    options.mHTL double = 1/500^1.5; % Suits simulations with only generalists
    options.bQuadraticHTL logical = false;
    options.bDecliningHTL logical = false;
end

L = 30;
T = 10;
n = length(range);
%
% Set up the parameter fields:
%
for i = 1:n
    pp(i) = p;
    if ~strcmp(sField, 'mortHTL') && ~strcmp(sField, 'sinking')
        eval( strcat('pp(i).',sField, '= range(i);') );
    end
end
%
% Do the simulations:
%
N = zeros(1,n);
Nlower = N;
Nupper = N;
B = zeros(n, p.nGroups);
Blower = B;
Bupper = B;
nGroups = p.nGroups;
mc = NaN;
Bc = NaN;


if options.bParallel
    parfor i = 1:n
        %
        % Set up simulation:
        %
        ppp = pp(i);
        ppp.tEnd = 1000;
        %
        % Simulate:
        %
        if strcmp(sField, 'mortHTL')
            setHTL(range(i), options.mHTL, options.bQuadraticHTL, options.bDecliningHTL);
        end

        if strcmp(sField, 'sinking')
            setSinkingPOM(ppp, range(i));
        end
        sim = simulateChemostat(ppp, L, T);
        %
        % Calc statistics:
        %
        ixAve = find(sim.t > sim.t(end)/2);

        N(i) = exp( mean( log(sim.N(ixAve))));
        Nlower(i) = min(sim.N(ixAve));
        Nupper(i) = max(sim.N(ixAve));

        if isfield(p,'idxSi')
            Si(i) = exp( mean( log(sim.Si(ixAve))));
            Silower(i) = min(sim.Si(ixAve));
            Siupper(i) = max(sim.Si(ixAve));
        end

        for iGroup = 1:nGroups
            ix = (ppp.ixStart(iGroup):ppp.ixEnd(iGroup)) - ppp.idxB + 1;
            B(i,iGroup) = max(1e-20, exp(mean(log(sum(sim.B(ixAve,ix),2)))));
            Blower(i,iGroup) = min(sum( sim.B(ixAve,ix),2 ));
            Bupper(i,iGroup) = max(sum( sim.B(ixAve,ix),2 ));
        end

        %Bpnm(i,:) = calcPicoNanoMicro(sim);
    end
else
    for i = 1:n
        runSim(i)
    end
end
%%
% Plot
%

%
% Nutrients:
%
if options.bPlotSpectra
    tiledlayout(1,2)
    nexttile
end

% Bifurcation plot of biomasses:
semilogy(range, N,'--','linewidth',2,'color',p.colNutrients{p.idxN});
hold on
patch([range range(end:-1:1)], [Nlower, Nupper(end:-1:1)], 0.75*[0,0,1], ...
    'edgecolor','none','facealpha',0.15)

legendentries(1) = semilogy(range, N,'--','linewidth',2,'color',p.colNutrients{p.idxN});
sLegend{1} = 'N';

if isfield(p,'idxSi')
    semilogy(range, Si,'--','linewidth',2,'color',p.colNutrients{p.idxN});
    patch([range range(end:-1:1)], [Silower, Siupper(end:-1:1)], 0.75*[0,0,1], ...
        'edgecolor','none','facealpha',0.15)

    legendentries(2) = semilogy(range, Si,'--','linewidth',2,'color',p.colNutrients{p.idxSi});
    sLegend{2} = 'Si';
end

%
% Biomass:
%
for iGroup = 1:nGroups
    patch([range range(end:-1:1)], [Blower(:,iGroup)', Bupper(end:-1:1,iGroup)'],...
        p.colGroup{iGroup}, 'edgecolor','none','facealpha',0.15)

    lwd = 1;
    if (p.typeGroups(iGroup) >= 10)
        lwd = max(1, 1 + log10( max(p.m(p.ixStart(iGroup):p.ixEnd(iGroup)))));
    end
    %% 
    legendentries(p.idxB-2+iGroup) = ...
        semilogy(range, B(:,iGroup),'color', p.colGroup{iGroup},...
        'linewidth',lwd);

    sLegend{iGroup+p.idxB-2} = p.nameGroup{iGroup};
end
set(gca,'xscale','log','yscale','log')
hold off
%
% Legend
%
legend(legendentries, sLegend, 'location','northwest','box','off')
xlabel(sField)
ylabel('Biomass ({\mu}g/l)')
ylim([1e-4 1e3])
xlim([min(range), max(range)])
%
% Spectra
%
if options.bPlotSpectra
    nexttile
    surface(mc,range,log10(Bc))
    shading flat
    set(gca,'xscale','log','yscale','log')
    ylabel(sField)
    xlabel('Cell mass ({\mu}g_C)')
    xlim('tight')
    ylim('tight')
    colorbar
end

    function runSim(i)
        %
        % Set up simulation:
        %
        ppp = pp(i);
        ppp.tEnd = 1000;
        %
        % Simulate:
        %
        if strcmp(sField, 'mortHTL')
            setHTL(range(i), options.mHTL, options.bQuadraticHTL, options.bDecliningHTL);
        end
        sim = simulateChemostat(ppp, L, T);
        %
        % Calc statistics:
        %
        ixAve = find(sim.t > sim.t(end)/2);

        N(i) = exp( mean( log(sim.N(ixAve))));
        Nlower(i) = min(sim.N(ixAve));
        Nupper(i) = max(sim.N(ixAve));

        if isfield(p,'idxSi')
            Si(i) = exp( mean( log(sim.Si(ixAve))));
            Silower(i) = min(sim.Si(ixAve));
            Siupper(i) = max(sim.Si(ixAve));
        end

        for iGroup = 1:nGroups
            ix = (ppp.ixStart(iGroup):ppp.ixEnd(iGroup)) - ppp.idxB + 1;
            B(i,iGroup) = exp(mean(log(sum(sim.B(ixAve,ix),2))));
            Blower(i,iGroup) = min(sum( sim.B(ixAve,ix),2 ));
            Bupper(i,iGroup) = max(sum( sim.B(ixAve,ix),2 ));
        end

        if isnan(mc)
            [mc, Bc] = calcCommunitySpectrum(sim.B, sim);
        else
            [mc, Bc(i,:)] = calcCommunitySpectrum(sim.B, sim);
        end
        %Bpnm(i,:) = calcPicoNanoMicro(sim);
    end
end