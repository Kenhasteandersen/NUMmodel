function sim=getSimRates(sim)

        sLibName = loadNUMmodelLibrary();
        ixTime = find(sim.t>(max(sim.t)-365)); %nTime = length(sim.t(sim.t >= max(sim.t-365))); % Just do the last year
        nZ = length(sim.z);
        % Find unicellular plankton groups:
        % [~,idxU]=find(sim.p.typeGroups<10);
        % ixUPlankton = sim.p.ixStart(idxU(1)):sim.p.ixEnd(idxU(end));
        
        for iTime = ixTime
            for k = 1:nZ
                if ~isnan(sim.N(iTime,k))
                    % Get the functions per volume at each depth and time:
                    if sim.p.nNutrients==3
                        u = [squeeze(sim.N(iTime,k)), ...
                            squeeze(sim.DOC(iTime,k)), ...
                            squeeze(sim.Si(iTime,k)), ...
                            squeeze(sim.B(iTime,k,:))'];
                    else
                        u = [squeeze(sim.N(iTime,k)), ...
                            squeeze(sim.DOC(iTime,k)), ...
                            squeeze(sim.B(iTime,k,:))'];
                    end
                    rates = getRates(sim.p, u, sim.L(iTime,k), sim.T(iTime,k), sLibName);
                    jLreal(iTime,k,:)=rates.jLreal;
                    jF(iTime,k,:)=rates.jF;
                    jDOC(iTime,k,:)=rates.jDOC;
                  
                end
            end
            sim.jLreal = jLreal;
            sim.jF = jF;
            sim.jDOC = jDOC;

        end