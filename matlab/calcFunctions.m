%
% Calculate the functions from any simulation. Only the last year is
% calculated for watercolumn and global simulations.
%
% In:
%   sim - simulation structure
%   bPrintSummary - prints out a summary
%
% Out:
%   sim - same as input but with fields of global function added:
%         sim.ProdGross - gross primary production (mgC/m2/d)
%         sim.ProdNet   - net primary production (mgC/m2/d)
%         sim.ProdHTL   - amount of carbon extracted from the HTL mortality
%         sim.Bpico, sim.Bnano, sim.Bmicro - pico, micro, and nano plankton
%                         biomasses (mgC/m2)
%         sim.Bplankton - All planktion, except (dead) POM (mgC/m2)
%         sim.ChlArea (mgChl/m2)
%         
%         For global simulations additional fields are:
%         sim.Ntotal    - total dissolved N as a function of time (mugN)
%         sim.DOCtotal  - total dissolved DOC (mugC)
%         sim.Btotal    - total biomass (mugC)
%         sim.ProdNetTotal - total NPP (mgC/d)
%
%         sim.ProdNetAnnual - annual average NPP (mgC/m2/d)
%         sim.ProdHTLAnnual - annual avearge HTL production (mgC/m2/d)
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
        sim.Btotal = sum(sim.B,2) * sim.p.widthProductiveLayer; % mgC/m^2
        sim.ChlArea = sum( calcChl( sim.B(end,:), sim.rates, sim.L) ) * sim.p.widthProductiveLayer;
        %ChlArea now in mgChl/m^2
        if options.bPrintSummary
            fprintf("----------------------------------------------\n")
            fprintf("Final total biomass:  %8.3f mgC/m2\n", sim.Btotal(end));
            fprintf("Final Chl:            %8.3f gChl/m2\n", sim.ChlArea)
            fprintf("Final gross PP:       %8.3f mgC/m2/day\n", sim.ProdGross)
            fprintf("Final net PP:         %8.3f mgC/m2/day\n", sim.ProdNet)
            fprintf("Final net bact prod.: %8.3f mgC/m2/day\n", sim.ProdBact)
            fprintf("Final HTL production: %8.3f mgC/m2/day\n", sim.ProdHTL)
            fprintf("----------------------------------------------\n")
        end

    case 'watercolumn'
        sLibName = loadNUMmodelLibrary();
        ixTime = find(sim.t>(max(sim.t)-365)); %nTime = length(sim.t(sim.t >= max(sim.t-365))); % Just do the last year
        nZ = length(sim.z);
        % Find plankton groups:
        if isfield(sim.p,'ixPOM')
            ixPlankton = sim.p.idxB:(sim.p.ixStart(sim.p.ixPOM)-1);
        else
            ixPlankton = sim.p.idxB:sim.p.n;
        end
        for iTime = ixTime
            ProdGross = 0;
            ProdNet = 0;
            ProdHTL = 0;
            ProdBact = 0;
            Bpico = 0;
            Bnano = 0;
            Bmicro = 0;
            Bplankton = 0;
            ChlArea = 0; 
            ChlVolume = zeros(1, nZ); 
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
                    [ProdGross1, ProdNet1,ProdHTL1,ProdBact1,~,Bpico1,Bnano1,Bmicro1] = ...
                        getFunctions(u, sim.L(iTime,k), sim.T(iTime,k), sLibName);
                    % Multiply by the thickness of each layer:
                    %rates = getRates(sim.p, u, sim.L(iTime,k), sim.T(iTime,k), sLibName);
                    
                    conv = sim.dznom(k);
                    ProdGross = ProdGross + ProdGross1*conv;
                    ProdNet = ProdNet + ProdNet1*conv;
                    ProdHTL = ProdHTL +ProdHTL1*conv;
                    ProdBact = ProdBact + ProdBact1*conv;
                    %eHTL = eHTL + eHTL1/length(sim.z);
                    Bpico = Bpico + Bpico1*conv; % mgC/m2
                    Bnano = Bnano + Bnano1*conv;
                    Bmicro = Bmicro + Bmicro1*conv;
                    Bplankton = Bplankton + sum(u(ixPlankton))*conv;
                    % Chl:
                    rates = getRates(sim.p, u, sim.L(iTime,k), sim.T(iTime,k), sLibName);
                    tmp =  calcChl( squeeze(sim.B(iTime,k,:)), rates, sim.L(iTime,k));
                    ChlArea = ChlArea + tmp .* sim.dznom(k); %tmp is already the sum over all sizes, so no need to do sum(tmp) again
                    ChlVolume(k) = tmp;
                end
            end

            iTimenow = iTime - ixTime(1)+1;
            sim.ProdGross(iTimenow) = ProdGross;
            sim.ProdNet(iTimenow) = ProdNet;
            sim.ProdHTL(iTimenow) = ProdHTL;
            sim.ProdBact(iTimenow) = ProdBact;
            %%sim.eHTL(i,j,iTime) = eHTL;
            sim.Bpico(iTimenow) = Bpico;
            sim.Bnano(iTimenow) = Bnano;
            sim.Bmicro(iTimenow) = Bmicro;
            sim.Bplankton(iTimenow) = Bplankton;
            sim.Btotal = sim.Bpico + sim.Bnano + sim.Bmicro; %mgC/m^2
            sim.ChlArea(iTimenow) = ChlArea; 
            sim.ChlVolume(iTimenow,:) = ChlVolume;
        end
        sim.eHTL = sim.ProdHTL./sim.ProdNet;
        sim.ePP = sim.ProdNet./sim.ProdGross;

        if options.bPrintSummary
            fprintf("----------------------------------------------\n")
            fprintf("Average plankton biomass: %9.3f mgC/m2\n", mean(sim.Bplankton));
            fprintf("Average Chl:              %9.3f gChl/m2\n", mean(sim.ChlArea)) %units corrected (there was a time dimension)
            fprintf("Average gross PP:         %9.3f mgC/m2/day\n", mean(sim.ProdGross))
            fprintf("Average net PP:           %9.3f mgC/m2/day\n", mean(sim.ProdNet))
            fprintf("Average HTL production:   %9.3f mgC/m2/day\n", mean(sim.ProdHTL))
            fprintf("----------------------------------------------\n")
        end

    case 'global'
        
        if ~isfield(sim,'ProdGross')
            sLibName = loadNUMmodelLibrary();
            ixTime = find(sim.t>(max(sim.t)-365)); %nTime = length(sim.t(sim.t >= max(sim.t-365))); % Just do the last year
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
            
            for iTime = ixTime
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
                                %Bphytotmp = Bphytotmp + calcBphyto()
                                % Chl:
                                rates = getRates(p, u, L(iTime,i,j,k), T(iTime,i,j,k), sLibName);
                                tmp =  calcChl( squeeze(B(iTime,i,j,k,:)), rates, L(iTime,i,j,k)) ;%/ 1000; % Convert to mg
                                if ~isnan(tmp)
                                    %ChlArea(i,j) = ChlArea(i,j) + tmp * dz(i,j,k);
                                    ChlVolume(iTime,i,j,k) = ChlVolume(iTime,i,j,k) + tmp; %mgChl/m^3
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
            % Global totals over the last year
            %
            calcTotal = @(u) sum(u(ix(:)).*dv(ix(:))); % mug/day

            for i = 1:length(ixTime)
                sim.Ntotal(i) = calcTotal(sim.N(i,:,:,:));
                sim.DOCtotal(i) = calcTotal(sim.DOC(i,:,:,:)); % mugC
                sim.Btotal(i) = 0;
                for j = 1:(sim.p.n-sim.p.idxB+1)
                    sim.Btotal(i) = sim.Btotal(i) + calcTotal(sim.B(i,:,:,:,j)); % mugC
                end

                sim.ProdNetTotal(i) = sum(sum(squeeze(sim.ProdNet(i,:,:)).*dx.*dy)); % mgC/day
                sim.ProdGrossTotal(i) = sum(sum(squeeze(sim.ProdGross(i,:,:)).*dx.*dy)); % mgC/day
                sim.ProdHTLTotal(i) = sum(sum(squeeze(sim.ProdHTL(i,:,:)).*dx.*dy));% mgC/day
            end

            sim.ChlTotal = ChlArea.*dx.*dy;
            sim.ChlTotal = sum(sim.ChlTotal(:));% mgChl
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
        % Annual global means:
        %
        if (~isfield(sim, 'ProdNetAnnual'))
            sim.ProdNetAnnual = mean(sim.ProdNet(ixTime,:,:),1);
            sim.ProdHTLAnnual = mean(sim.ProdHTL(ixTime,:,:),1);
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

function Bphyto = calcBphyto(B, rates) % Calculate the phytoplankton biomass
    Bphyto = sum( B .* rates.jLreal./(rates.jLreal+rates.jFreal+rates.JDOCreal) );
end

end