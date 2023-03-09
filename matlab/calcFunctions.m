%
% Calculate the functions from any simulation.
%
% In:
%   sim - simulation structure
%   bPrintSummary - prints out a summary
%
% Out:
%   sim - same as input but with fields of global function added:
%         sim.ProdGross - gross primary production (gC/m2/yr)
%         sim.ProdNet   - net primary production (gC/m2/yr)
%         sim.ProdHTL   - amount of carbon extracted from the HTL mortality
%         sim.Bpico, sim.Bnano, sim.Bmicro - pico, micro, and nano plankton
%                         biomasses (gC/m2)
%         sim.ChlArea (gChl/m2)
%
%         For global simulations additional fields are:
%         sim.Ntotal    - total dissolved N as a function of time (mugN)
%         sim.DOCtotal  - total dissolved DOC (mugC)
%         sim.Btotal    - total biomass (mugC)
%         sim.ProdNetTotal - total NPP (gC/yr)
%
%         sim.ProdNetAnnual - annual average NPP (gC/m2/yr)
%         sim.ProdHTLAnnual - annual avearge HTL production (gc/m2/yr)
%
function sim = calcFunctions(sim, options)

arguments
    sim struct;
    options.bPrintSummary = true;
end


switch sim.p.nameModel

    case 'chemostat'
        if sim.p.nNutrients==3
            u = [sim.N(end), sim.DOC(end),sim.Si(end), sim.B(end,:)];
            [sim.ProdGross, sim.ProdNet, sim.ProdHTL, sim.ProdBact, sim.eHTL,...
                sim.Bpico, sim.Bnano, sim.Bmicro] = ...=
                getFunctions(u, sim.L, sim.T);
        else
            u = [sim.N(end), sim.DOC(end), sim.B(end,:)];
            [sim.ProdGross, sim.ProdNet, sim.ProdHTL, sim.ProdBact, sim.eHTL,...
                sim.Bpico, sim.Bnano, sim.Bmicro] = ...=
                getFunctions(u, sim.L, sim.T);
        end
        % Multiply by the assumed depth of the productive layer:
        sim.ProdGross = sim.ProdGross * sim.p.widthProductiveLayer;
        sim.ProdNet = sim.ProdNet * sim.p.widthProductiveLayer;
        sim.ProdHTL = sim.ProdHTL * sim.p.widthProductiveLayer;
        sim.ProdBact = sim.ProdBact * sim.p.widthProductiveLayer;
        sim.Bpico = sim.Bpico * sim.p.widthProductiveLayer;
        sim.Bnano = sim.Bnano * sim.p.widthProductiveLayer;
        sim.Bmicro = sim.Bmicro * sim.p.widthProductiveLayer;
        sim.Btot = sum(sim.B,2) * sim.p.widthProductiveLayer/1000;
        sim.ChlArea = sum( calcChl( sim.B(end,:), sim.rates, sim.L) ) * sim.p.widthProductiveLayer/ 1000; % Convert to g

        if options.bPrintSummary
            fprintf("----------------------------------------------\n")
            fprintf("Final total biomass:  %.2f gC/m2\n", sim.Btot(end));
            fprintf("Final Chl:            %.2f gChl/m2\n", sim.ChlArea)
            fprintf("Final gross PP:       %.2f gC/m2/yr\n", sim.ProdGross)
            fprintf("Final net PP:         %.2f gC/m2/yr\n", sim.ProdNet)
            fprintf("Final HTL production: %.2f gC/m2/yr\n", sim.ProdHTL)
            fprintf("----------------------------------------------\n")
        end

    case 'watercolumn'
        nTime = length(sim.t);
        nZ = length(sim.z);
        ChlArea = zeros(nTime,1);
        ChlVolume = zeros(nTime, nZ);
        jLreal = zeros(nTime,nZ,sim.p.n-sim.p.idxB+1);
        for iTime = 1:nTime
            ProdGross = 0;
            ProdNet = 0;
            ProdHTL = 0;
            ProdBact = 0;
            Bpico = 0;
            Bnano = 0;
            Bmicro = 0;
            % Integrate over depth:
            for k = 1:nZ
                if ~isnan(sim.N(k,iTime))
                    % Get the functions per volume at each depth and time:
                    if sim.p.nNutrients==3
                        u = [squeeze(sim.N(k,iTime)), ...
                            squeeze(sim.DOC(k,iTime)), ...
                            squeeze(sim.Si(k,iTime)), ...
                            squeeze(sim.B(k,:,iTime))];
                    else
                        u = [squeeze(sim.N(k,iTime)), ...
                            squeeze(sim.DOC(k,iTime)), ...
                            squeeze(sim.B(k,:,iTime))];
                    end
                    [ProdGross1, ProdNet1,ProdHTL1,ProdBact1,~,Bpico1,Bnano1,Bmicro1] = ...
                        getFunctions(u, sim.L(k,iTime), sim.T(k,iTime));
                    % Multiply by the thickness of each layer:
                    conv = sim.dznom(k);
                    ProdGross = ProdGross + ProdGross1*conv;
                    ProdNet = ProdNet + ProdNet1*conv;
                    ProdHTL = ProdHTL +ProdHTL1*conv;
                    ProdBact = ProdBact + ProdBact1*conv;
                    %eHTL = eHTL + eHTL1/length(sim.z);
                    Bpico = Bpico + Bpico1*conv; % gC/m2
                    Bnano = Bnano + Bnano1*conv;
                    Bmicro = Bmicro + Bmicro1*conv;
                    % Chl:
                    rates = getRates(sim.p, u, sim.L(k,iTime), sim.T(k,iTime));
                    tmp =  calcChl( squeeze(sim.B(k,:,iTime)), rates, sim.L(k,iTime)) / 1000; % Convert to g
                    if ~isnan(tmp)
                        jLreal(iTime,k,:) = rates.jLreal;
                        ChlArea(iTime) = ChlArea(iTime) + sum(tmp) * sim.dznom(k);
                        ChlVolume(iTime,k) = ChlVolume(iTime,k) + sum(tmp);
                    end
                end
            end
            sim.ProdGross(iTime) = ProdGross;
            sim.ProdNet(iTime) = ProdNet;
            sim.ProdHTL(iTime) = ProdHTL;
            sim.ProdBact(iTime) = ProdBact;
            %%sim.eHTL(i,j,iTime) = eHTL;
            sim.Bpico(iTime) = Bpico;
            sim.Bnano(iTime) = Bnano;
            sim.Bmicro(iTime) = Bmicro;
        end
        sim.eHTL = sim.ProdHTL./sim.ProdNet;
        sim.ePP = sim.ProdNet./sim.ProdGross;
        sim.ChlArea = ChlArea;
        sim.ChlVolume = ChlVolume;
        sim.jLreal = jLreal;

        if options.bPrintSummary
            fprintf("----------------------------------------------\n")
            fprintf("Average total biomass:  %.2f gC/m2\n", mean(sim.Bpico+sim.Bnano+sim.Bmicro));
            fprintf("Average Chl:            %.2f gChl/m2/yr\n", mean(sim.ChlArea))
            fprintf("Average gross PP:       %.2f gC/m2/yr\n", mean(sim.ProdGross))
            fprintf("Average net PP:         %.2f gC/m2/yr\n", mean(sim.ProdNet))
            fprintf("Average HTL production: %.2f gC/m2/yr\n", mean(sim.ProdHTL))
            fprintf("----------------------------------------------\n")
        end

    case 'global'
        %
        % Primary production in gC per m2/year::
        %
        if ~isfield(sim,'ProdGross')
            % Get grid volumes:
            load(sim.p.pathGrid,'dv','dz','dx','dy');
            ix = ~isnan(sim.N(:,:,1,1)); % Find all relevant grid cells

            sim.ProdGross = zeros(length(sim.x), length(sim.y), length(sim.t));
            sim.ProdNet = sim.ProdGross;
            sim.ProdHTL = sim.ProdGross;
            sim.ProdBact = sim.ProdGross;
            sim.Bpico = sim.ProdGross;
            sim.Bnano = sim.ProdGross;
            sim.Bmicro = sim.ProdGross;
            ChlArea = 0*dx;
            ChlVolume = zeros(length(sim.t),length(sim.x), length(sim.y), length(sim.z));

            nTime = length(sim.t);
            nX = length(sim.x);
            nY = length(sim.y);
            nZ = length(sim.z);
            for iTime = 1:nTime
                for i = 1:nX
                    for j = 1:nY
                        ProdGross = 0;
                        ProdNet = 0;
                        ProdHTL = 0;
                        ProdBact = 0;
                        Bpico = 0;
                        Bnano = 0;
                        Bmicro = 0;
                        for k = 1:nZ
                            if ~isnan(sim.N(i,j,k,iTime))
                                if isfield(sim.p,'idxSi')
                                    u = [squeeze(sim.N(i,j,k,iTime)), ...
                                        squeeze(sim.DOC(i,j,k,iTime)), ...
                                        squeeze(sim.Si(i,j,k,iTime)), ...
                                        squeeze(sim.B(i,j,k,:,iTime))'];
                                else
                                    u = [squeeze(sim.N(i,j,k,iTime)), ...
                                        squeeze(sim.DOC(i,j,k,iTime)), ...
                                        squeeze(sim.B(i,j,k,:,iTime))'];
                                end
                                [ProdGross1, ProdNet1,ProdHTL1,ProdBact1, ~,Bpico1,Bnano1,Bmicro1] = ...
                                    getFunctions(u, sim.L(i,j,k,iTime), sim.T(i,j,k,iTime));
                                conv = squeeze(dz(i,j,k));
                                ProdGross = ProdGross + ProdGross1*conv; % gC/m2/yr
                                ProdNet = ProdNet + ProdNet1*conv;
                                ProdHTL = ProdHTL +ProdHTL1*conv;
                                ProdBact = ProdBact +ProdBact1*conv;
                                %eHTL = eHTL + eHTL1/length(sim.z);
                                Bpico = Bpico + Bpico1*dz(i,j,k); % gC/m2
                                Bnano = Bnano + Bnano1*dz(i,j,k);
                                Bmicro = Bmicro + Bmicro1*dz(i,j,k);
                                % Chl:
                                rates = getRates(sim.p, u, sim.L(i,j,k,iTime), sim.T(i,j,k,iTime));
                                tmp =  calcChl( squeeze(sim.B(i,j,k,:,iTime)), rates, sim.L(i,j,k,iTime)) ;%/ 1000; % Convert to mg
                                if ~isnan(tmp)
                                    ChlArea(i,j) = ChlArea(i,j) + tmp * dz(i,j,k);
                                    ChlVolume(iTime,i,j,k) = ChlVolume(iTime,i,j,k) + tmp;
                                end
                            end
                        end
                        sim.ProdGross(i,j,iTime) = ProdGross;
                        sim.ProdNet(i,j,iTime) = ProdNet;
                        sim.ProdHTL(i,j,iTime) = ProdHTL;
                        sim.ProdBact(i,j,iTime) = ProdBact;
                        sim.Bpico(i,j,iTime) = Bpico;
                        sim.Bnano(i,j,iTime) = Bnano;
                        sim.Bmicro(i,j,iTime) = Bmicro;
                    end
                end
            end
            sim.eHTL = sim.ProdHTL./sim.ProdNet;
            sim.ePP = sim.ProdNet./sim.ProdGross;
            sim.ChlArea = ChlArea/length(sim.t);
            sim.ChlVolume = ChlVolume/length(sim.t);
            %
            % Global totals
            %
            calcTotal = @(u) sum(u(ix(:)).*dv(ix(:))); % mug/day

            for i = 1:length(sim.t)
                sim.Ntotal(i) = calcTotal(sim.N(:,:,:,i));
                sim.DOCtotal(i) = calcTotal(sim.DOC(:,:,:,i)); % mugC
                sim.Btotal(i) = 0;
                for j = 1:(sim.p.n-sim.p.idxB+1)
                    sim.Btotal(i) = sim.Btotal(i) + calcTotal(sim.B(:,:,:,j,i)); % mugC
                end

                sim.ProdNetTotal(i) = sum(sum(sim.ProdNet(:,:,i).*dx.*dy)); % gC/yr
                sim.ProdGrossTotal(i) = sum(sum(sim.ProdGross(:,:,i).*dx.*dy)); % gC/yr
                sim.ProdHTLTotal(i) = sum(sum(sim.ProdHTL(:,:,i).*dx.*dy));
            end

            sim.ChlTotal = ChlArea.*dx.*dy;
            sim.ChlTotal = sum(sim.ChlTotal(:));%/ 1000; %gChl
        end
        %
        % Watercolumn totals:
        %
        % for i = 1:length(sim.x)
        %     for j = 1:length(sim.y)
        %         for k = 1:length(sim.t)
        %             % gC per m2/year:
        %             sim.ProdNetPerArea(i,j,k) = sum(squeeze(sim.ProdNet(i,j,:,k)).*squeeze(dz(i,j,:)))*365*10*10*1000;
        %         end
        %     end
        % end
        %
        % Annual global totals. Less accurate than calculating them via the flag
        % "bCalcGlobalAnnual" in simulateGlobal()
        %
        if (~isfield(sim, 'ProdNetAnnual'))
            sim.ProdNetAnnual = mean(sim.ProdNet,3);
            sim.ProdHTLAnnual(:,:,i) = mean(sim.ProdHTL,3);
            %zeros(length(sim.x), length(sim.y), floor(sim.t(end)/365));
            %for i = 1:sim.t(end)/365
            %    ixTime = sim.t>365*(i-1) & sim.t<=365*i;
            %    sim.ProdNetTotalAnnual = mean(sim.ProdNetTotal)/length(ixTime); % mugC/l/day
            %    sim.ProdNetAnnual(:,:,i) = /length(ixTime); % gC/m2/yr
            %    length(ixTime); % gC/m2/yr
            %end
            %sim.ProdNetTotalAnnual = sim.ProdNetTotalAnnual/1000/1000/1000/1000; % GTon carbon/year
        end

    otherwise
        disp('Model unknown; functions not calculated');

end
sim.Btotal = sim.Bpico + sim.Bnano + sim.Bmicro;
