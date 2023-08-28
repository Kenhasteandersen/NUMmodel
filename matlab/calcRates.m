%
% Calculate jLreal, jFreal & jDOC for a simulation
%
% If tDay<0, average over the last month of the simulation
%
% In:
%    sim - simulation structure
%    tDay - time (-1 by default)
%
% Out :
%   Rates - jLreal, jFreal & jDOC 
%

function Rates = calcRates(sim,tDay)

arguments
    sim struct
    tDay double = sim.t(end);
end

if tDay<0
    %Average over the last month
    if sim.t(end)<30 %if the simulation lasted less than one month then the average is done over all the simulation
        time = 1;
    else
        [~, time] = min(abs(sim.t-(sim.t(end)-29)));
    end
    time = time:length(sim.t);
else
    [~, time] = min(abs(sim.t-tDay));
end

%
%Initialisation
%
switch sim.p.nameModel
    case 'chemostat'
        jLreal = zeros(length(time),sim.p.n-sim.p.nNutrients);
    case 'watercolumn'
        jLreal = zeros(length(time),length(sim.z),sim.p.n-sim.p.nNutrients);
    case 'global'
        jLreal = zeros(length(time),length(sim.x),length(sim.y),length(sim.z),sim.p.n-sim.p.nNutrients);
end
jFreal = jLreal; jDOC = jLreal;

%
% Calculation
%
for i = 1:length(time)

    iTime = time(i);

    switch sim.p.nameModel

        case 'chemostat'
            rates = getRates(sim.p, sim.u(iTime,:), mean(sim.L), sim.T);
            jLreal(i,:) = rates.jLreal;
            jFreal(i,:) = rates.jFreal;
            jDOC(i,:) = rates.jDOC;
        
        case 'watercolumn'
            for iDepth = 1:length(sim.z)
                if isfield(sim,'Si')
                    u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), sim.Si(iTime,iDepth), ...
                        squeeze(sim.B(iTime,iDepth,:))'];
                else
                    u = [sim.N(iTime,iDepth), sim.DOC(iTime,iDepth), squeeze(sim.B(iTime,iDepth,:))'];
                end
                rates = getRates(sim.p, u, sim.L(iTime,iDepth), sim.T(iTime,iDepth));
                jLreal(i,iDepth,:) = rates.jLreal;
                jFreal(i,iDepth,:) = rates.jFreal;
                jDOC(i,iDepth,:) = rates.jDOC;
            end
    
        case 'global'
            for ixX = 1:length(sim.x)
                for ixY = 1:length(sim.y)
                    for iDepth = 1:length(sim.z)
                        if isfield(sim,'Si')
                            u = [sim.N(iTime,ixX, ixY,iDepth), ...
                                sim.DOC(iTime,ixX, ixY,iDepth), ...
                                sim.Si(iTime,ixX, ixY,iDepth), ...
                                squeeze(sim.B(iTime,ixX, ixY, iDepth, :))'];
                        else
                            u = [sim.N(iTime,ixX, ixY,iDepth), ...
                                sim.DOC(iTime,ixX, ixY,iDepth), ...
                                squeeze(sim.B(iTime,ixX, ixY, iDepth, :))'];
                        end
                        rates = getRates(sim.p, u, sim.L(iTime,ixX,ixY,iDepth), sim.T(iTime,ixX,ixY,iDepth));
                        jLreal(i,ixX,ixY,iDepth,:) = rates.jLreal;
                        jFreal(i,ixX,ixY,iDepth,:) = rates.jFreal;
                        jDOC(i,ixX,ixY,iDepth,:) = rates.jDOC;
                    end
                end
            end
    end
end

Rates.jLreal = squeeze(mean(jLreal,1));
Rates.jFreal = squeeze(mean(jFreal,1));
Rates.jDOC = squeeze(mean(jDOC,1));

