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
%         sim.ProdNet   - net primary production (mgC/m2/d)
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
            fprintf("Final total biomass:  %8.3f gC/m2\n", sim.Btot(end));
            fprintf("Final Chl:            %8.3f gChl/m2\n", sim.ChlArea)
            fprintf("Final gross PP:       %8.3f gC/m2/yr\n", sim.ProdGross)
            fprintf("Final net PP:         %8.3f gC/m2/yr\n", sim.ProdNet)
            fprintf("Final net bact prod.: %8.3f gC/m2/yr\n", sim.ProdBact)
            fprintf("Final HTL production: %8.3f gC/m2/yr\n", sim.ProdHTL)
            fprintf("----------------------------------------------\n")
        end

    case 'watercolumn'
        sLibName = loadNUMmodelLibrary();
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
            Nitr = zeros(nTime);
            % Integrate over depth:
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
                   
                    if k<=9
                        Nitr(iTime) = Nitr(iTime) + u(1) * sim.dznom(k);
                    end


                    [ProdGross1, ProdNet1,ProdHTL1,ProdBact1,~,Bpico1,Bnano1,Bmicro1] = ...
                        getFunctions(u, sim.L(iTime,k), sim.T(iTime,k), sLibName);
                    % Multiply by the thickness of each layer:
                    conv = sim.dznom(k);
                    ProdGross = ProdGross + ProdGross1*conv;
                    ProdNet = ProdNet + ProdNet1*conv; % mgC/m2/d
                    ProdHTL = ProdHTL +ProdHTL1*conv;
                    ProdBact = ProdBact + ProdBact1*conv;
                    %eHTL = eHTL + eHTL1/length(sim.z);
                    Bpico = Bpico + Bpico1*conv; % gC/m2
                    Bnano = Bnano + Bnano1*conv;
                    Bmicro = Bmicro + Bmicro1*conv;
                    % Chl:
                    rates = getRates(sim.p, u, sim.L(iTime,k), sim.T(iTime,k), sLibName);
                    tmp =  calcChl( squeeze(sim.B(iTime,k,:)), rates, sim.L(iTime,k)) / 1000; % Convert to g
                    if ~isnan(tmp)
                        jLreal(iTime,k,:) = rates.jLreal;
                        ChlArea(iTime) = ChlArea(iTime) + sum(tmp) * sim.dznom(k);
                        ChlVolume(iTime,k) = ChlVolume(iTime,k) + sum(tmp);
                    end
                    ratePOM(iTime, k, :) = rates.jPOM;
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
        sim.ratePOM = ratePOM;

        if options.bPrintSummary
            fprintf("----------------------------------------------\n")
            fprintf("Average total biomass:  %8.3f gC/m2\n", mean(sim.Bpico+sim.Bnano+sim.Bmicro));
            fprintf("Average Chl:            %8.3f gChl/m2/yr\n", mean(sim.ChlArea))
            fprintf("Average gross PP:       %8.3f gC/m2/yr\n", mean(sim.ProdGross))
            fprintf("Average net PP:         %8.3f mgC/m2/d\n", mean(sim.ProdNet)) % changed-- the others should too
            fprintf("Average HTL production: %8.3f gC/m2/yr\n", mean(sim.ProdHTL))
            fprintf("----------------------------------------------\n")
        end

    case 'global'
        
        if ~isfield(sim,'ProdGross')
            sLibName = loadNUMmodelLibrary();
            % Get grid volumes:
            load(sim.p.pathGrid,'dv','dz','dx','dy');
            ix = ~isnan(sim.N(1,:,:,1)); % Find all relevant grid cells

            ProdGross = zeros(length(sim.t), length(sim.x), length(sim.y));
            ProdNet = ProdGross;
            ProdHTL = ProdGross;
            ProdBact = ProdGross;
            Bpico = ProdGross;
            Bnano = ProdGross;
            Bmicro = ProdGross;
            ChlArea = 0*dx;
            ChlVolume = zeros(length(sim.t),length(sim.x), length(sim.y), length(sim.z));

            nTime = length(sim.t);
            nX = length(sim.x);
            nY = length(sim.y);
            nZ = length(sim.z);
            %
            % Extract fields from sim:
            %
            N = sim.N;
            DOC = sim.DOC;
            if isfield(sim.p,'idxSi')
                Si = sim.Si;
            else
                Si = 0;
            end
            B = sim.B;
            L = sim.L;
            T = sim.T;
            p = sim.p;
            
            for iTime = 1:nTime
                for i = 1:nX
                    for j = 1:nY
                        ProdGrosstmp = 0;
                        ProdNettmp = 0;
                        ProdHTLtmp = 0;
                        ProdBacttmp = 0;
                        Bpicotmp = 0;
                        Bnanotmp = 0;
                        Bmicrotmp = 0;
                        for k = 1:nZ
                            if ~isnan(N(iTime,i,j,k))
                                if isfield(p,'idxSi')
                                    u = [squeeze(N(iTime,i,j,k)), ...
                                        squeeze(DOC(iTime,i,j,k)), ...
                                        squeeze(Si(iTime,i,j,k)), ...
                                        squeeze(B(iTime,i,j,k,:))'];
                                else
                                    u = [squeeze(N(iTime,i,j,k)), ...
                                        squeeze(DOC(iTime,i,j,k)), ...
                                        squeeze(B(iTime,i,j,k,:))'];
                                end
                                [ProdGross1, ProdNet1,ProdHTL1,ProdBact1, ~,Bpico1,Bnano1,Bmicro1] = ...
                                    getFunctions(u, L(iTime,i,j,k), T(iTime,i,j,k), sLibName);
                                conv = squeeze(dz(i,j,k));
                                ProdGrosstmp = ProdGrosstmp + ProdGross1*conv; % gC/m2/yr
                                ProdNettmp = ProdNettmp + ProdNet1*conv;
                                ProdHTLtmp = ProdHTLtmp +ProdHTL1*conv;
                                ProdBacttmp = ProdBacttmp +ProdBact1*conv;
                                %eHTL = eHTL + eHTL1/length(sim.z);
                                Bpicotmp = Bpicotmp + Bpico1*dz(i,j,k); % gC/m2
                                Bnanotmp = Bnanotmp + Bnano1*dz(i,j,k);
                                Bmicrotmp = Bmicrotmp + Bmicro1*dz(i,j,k);
                                % Chl:
                                rates = getRates(p, u, L(iTime,i,j,k), T(iTime,i,j,k), sLibName);
                                tmp =  calcChl( squeeze(B(iTime,i,j,k,:)), rates, L(iTime,i,j,k)) ;%/ 1000; % Convert to mg
                                if ~isnan(tmp)
                                    %ChlArea(i,j) = ChlArea(i,j) + tmp * dz(i,j,k);
                                    ChlVolume(iTime,i,j,k) = ChlVolume(iTime,i,j,k) + tmp;
                                end
                            end
                        end
                        ProdGross(iTime,i,j) = ProdGrosstmp;
                        ProdNet(iTime,i,j) = ProdNettmp;
                        ProdHTL(iTime,i,j) = ProdHTLtmp;
                        ProdBact(iTime,i,j) = ProdBacttmp;
                        Bpico(iTime,i,j) = Bpicotmp;
                        Bnano(iTime,i,j) = Bnanotmp;
                        Bmicro(iTime,i,j) = Bmicrotmp;
                    end
                end
            end
            sim.ProdGross = ProdGross;
            sim.ProdNet = ProdNet;
            sim.ProdHTL = ProdHTL;
            sim.ProdBact = ProdBact;
            sim.Bpico = Bpico;
            sim.Bnano = Bnano;
            sim.Bmicro = Bmicro;
            sim.eHTL = sim.ProdHTL./sim.ProdNet;
            sim.ePP = sim.ProdNet./sim.ProdGross;
            sim.ChlArea = ChlArea/length(sim.t);
            sim.ChlVolume = ChlVolume/length(sim.t);
            %
            % Global totals
            %
            calcTotal = @(u) sum(u(ix(:)).*dv(ix(:))); % mug/day

            for i = 1:length(sim.t)
                sim.Ntotal(i) = calcTotal(sim.N(i,:,:,:));
                sim.DOCtotal(i) = calcTotal(sim.DOC(i,:,:,:)); % mugC
                sim.Btotal(i) = 0;
                for j = 1:(sim.p.n-sim.p.idxB+1)
                    sim.Btotal(i) = sim.Btotal(i) + calcTotal(sim.B(i,:,:,:,j)); % mugC
                end

                sim.ProdNetTotal(i) = sum(sum(squeeze(sim.ProdNet(i,:,:)).*dx.*dy)); % gC/yr
                sim.ProdGrossTotal(i) = sum(sum(squeeze(sim.ProdGross(i,:,:)).*dx.*dy)); % gC/yr
                sim.ProdHTLTotal(i) = sum(sum(squeeze(sim.ProdHTL(i,:,:)).*dx.*dy));
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
            sim.ProdNetAnnual = mean(sim.ProdNet,1);
            sim.ProdHTLAnnual(i,:,:) = mean(sim.ProdHTL,1);
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