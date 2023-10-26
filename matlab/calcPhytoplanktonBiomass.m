function Bphyto=calcPhytoplanktonBiomass(sim)

        sLibName = loadNUMmodelLibrary();
        ixTime = find(sim.t>(max(sim.t)-365)); %nTime = length(sim.t(sim.t >= max(sim.t-365))); % Just do the last year
        nZ = length(sim.z);
        % Find unicellular plankton groups:
        [~,idxU]=find(sim.p.typeGroups<10);
        ixUPlankton = sim.p.ixStart(idxU(1)):sim.p.ixEnd(idxU(end));
        
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
                    % Photosynthetic biomass
                    Bphyto(iTime,k,:)=calcBphyto(squeeze(sim.B(iTime,k,:)), rates,ixUPlankton);
                    
                    % tmp=findPhyto_ix(rates,ixUPlankton);
                    % ix_phyto(iTime,k,:)=findPhyto_ix(rates,ixUPlankton);
                end
            end
            iTimenow = iTime - ixTime(1)+1;
            % sim.ixPhyto(iTimenow,:,:) = ix_phyto;
            % sim.Bphyto=Bphyto;
        end
% Bpnmphyto = [Bpico, Bnano, Bmicro];

% I must set a ratio of what is autotroph
function Bphyto = calcBphyto(B,rates,ix) %check units!
        Bphyto_frac = rates.jLreal(ix)./(rates.jLreal(ix)+rates.jF(ix)+rates.jDOC(ix));
        Bphyto = rates.jLreal(ix)./(rates.jLreal(ix)+rates.jF(ix)+rates.jDOC(ix)).*B(ix);
        if isnan(Bphyto)
            Bhyto=0;
        end
        % ixPhyto=find(Bphyto_frac>=0.6);
        % Bphyto = rates.jLreal(ixPhyto)./(rates.jLreal(ixPhyto)+rates.jF(ixPhyto)+rates.jDOC(ixPhyto)).*B(ixPhyto);
        % Bphyto_true = Bphyto_frac(ixPhyto).*B(ixPhyto);

% function ixPhyto = findPhyto_ix(rates,ix)
% 
% Bphyto_frac = rates.jLreal(ix)./(rates.jLreal(ix)+rates.jF(ix)+rates.jDOC(ix));
%         % Bphyto = rates.jLreal(ix)./(rates.jLreal(ix)+rates.jF(ix)+rates.jDOC(ix)).*B(ix);
%         ixPhyto=find(Bphyto_frac>=0.6);