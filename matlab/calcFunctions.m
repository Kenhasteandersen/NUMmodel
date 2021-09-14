%
% Calculate the functions from any simulation.
%
% In:
%   sim - simulation structure
%
% Out:
%   sim - same as input but with fields of global function added:
%         sim.ProdGross - gross primary production (gC/m2/yr)
%         sim.ProdNet   - net primary production (gC/m2/yr)
%         sim.ProdHTL   - amount of carbon extracted from the HTL mortality
%         sim.Bpico, sim.Bnano, sim.Bmicro - pico, micro, and nano plankton
%                         biomasses (gC/m2/yr)
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
function sim = calcFunctions(sim)

switch sim.p.nameModel
    
    case 'chemostat'
        u = [sim.N(end), sim.DOC(end), sim.B(end,:)];
        [sim.ProdGross, sim.ProdNet, sim.ProdHTL, sim.eHTL,...
            sim.Bpico, sim.Bnano, sim.Bmicro] = ...
            getFunctions(u, sim.L, sim.T);
        % Multiply by the assumed depth of the productive layer:
        sim.ProdGross = sim.ProdGross * sim.p.depthProductiveLayer;
        sim.ProdNet = sim.ProdNet * sim.p.depthProductiveLayer;
        sim.ProdHTL = sim.ProdHTL * sim.p.depthProductiveLayer;
        sim.Bpico = sim.Bpico * sim.p.depthProductiveLayer;
        sim.Bnano = sim.Bnano * sim.p.depthProductiveLayer;
        sim.Bmicro = sim.Bmicro * sim.p.depthProductiveLayer;
        
    case 'watercolumn'
        nTime = length(sim.t);
        nZ = length(sim.z);
        for iTime = 1:nTime
            ProdGross = 0;
            ProdNet = 0;
            ProdHTL = 0;
            Bpico = 0;
            Bnano = 0;
            Bmicro = 0;
            % Integrate over depth:
            for k = 1:nZ
                if ~isnan(sim.N(k,iTime))
                    % Get the functions per volume at each depth and time:
                    u = [squeeze(sim.N(k,iTime)), ...
                        squeeze(sim.DOC(k,iTime)), ...
                        squeeze(sim.B(k,:,iTime))];
                    [ProdGross1, ProdNet1,ProdHTL1,eHTL,Bpico1,Bnano1,Bmicro1] = ...
                        getFunctions(u, sim.L(k,iTime), sim.T(k,iTime));
                    % Multiply by the thickness of each layer:
                    conv = sim.dznom(k);
                    ProdGross = ProdGross + ProdGross1*conv;
                    ProdNet = ProdNet + ProdNet1*conv;
                    ProdHTL = ProdHTL +ProdHTL1*conv;
                    %eHTL = eHTL + eHTL1/length(sim.z);
                    Bpico = Bpico + Bpico1*conv; % gC/m2
                    Bnano = Bnano + Bnano1*conv;
                    Bmicro = Bmicro + Bmicro1*conv;
                end
            end
            sim.ProdGross(iTime) = ProdGross;
            sim.ProdNet(iTime) = ProdNet;
            sim.ProdHTL(iTime) = ProdHTL;
            %%sim.eHTL(i,j,iTime) = eHTL;
            sim.Bpico(iTime) = Bpico;
            sim.Bnano(iTime) = Bnano;
            sim.Bmicro(iTime) = Bmicro;
        end
        
        
    case 'global'
        % Get grid volumes:
        load(sim.p.pathGrid,'dv','dz');
        ix = ~isnan(sim.N(:,:,1,1)); % Find all relevant grid cells
        
        %
        % Primary production in gC per m2/year::
        %
        if ~isfield(sim,'ProdGross')
            sim.ProdGross = zeros(length(sim.x), length(sim.y), length(sim.t));
            sim.ProdNet = sim.ProdGross;
            sim.ProdHTL = sim.ProdGross;
            sim.Bpico = sim.ProdGross;
            sim.Bnano = sim.ProdGross;
            sim.Bmicro = sim.ProdGross;
            
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
                        Bpico = 0;
                        Bnano = 0;
                        Bmicro = 0;
                        for k = 1:nZ
                            if ~isnan(sim.N(i,j,k,iTime))
                                u = [squeeze(sim.N(i,j,k,iTime)), ...
                                    squeeze(sim.DOC(i,j,k,iTime)), ...
                                    squeeze(sim.B(i,j,k,:,iTime))'];
                                [ProdGross1, ProdNet1,ProdHTL1,eHTL,Bpico1,Bnano1,Bmicro1] = ...
                                    getFunctions(u, sim.L(i,j,k,iTime), sim.T(i,j,k,iTime));
                                conv = squeeze(dz(i,j,k));
                                ProdGross = ProdGross + ProdGross1*conv;
                                ProdNet = ProdNet + ProdNet1*conv;
                                ProdHTL = ProdHTL +ProdHTL1*conv;
                                %eHTL = eHTL + eHTL1/length(sim.z);
                                Bpico = Bpico + Bpico1*dz(i,j,k); % gC/m2
                                Bnano = Bnano + Bnano1*dz(i,j,k);
                                Bmicro = Bmicro + Bmicro1*dz(i,j,k);
                            end
                        end
                        sim.ProdGross(i,j,iTime) = ProdGross;
                        sim.ProdNet(i,j,iTime) = ProdNet;
                        sim.ProdHTL(i,j,iTime) = ProdHTL;
                        %%sim.eHTL(i,j,iTime) = eHTL;
                        sim.Bpico(i,j,iTime) = Bpico;
                        sim.Bnano(i,j,iTime) = Bnano;
                        sim.Bmicro(i,j,iTime) = Bmicro;
                    end
                end
            end
        end
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
            sim.ProdNetTotal(i) = calcTotal(sim.ProdNet(:,:,i)); % gC/yr
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
