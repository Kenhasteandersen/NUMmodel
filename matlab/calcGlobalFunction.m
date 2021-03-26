function sim = calcGlobalFunction(sim)

loadNUMmodelLibrary();
calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(sim.p.mAdult)), sim.p.mAdult)
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
    
    for iTime = 1:length(sim.t)
        for i = 1:length(sim.x)
            for j = 1:length(sim.y)
                ProdGross = 0;
                ProdNet = 0;
                ProdHTL = 0;
                Bpico = 0;
                Bnano = 0;
                Bmicro = 0;
                for k = 1:length(sim.z)
                    if ~isnan(sim.N(i,j,k,iTime))
                        u = [squeeze(sim.N(i,j,k,iTime)), ...
                            squeeze(sim.DOC(i,j,k,iTime)), ...
                            squeeze(sim.B(i,j,k,:,iTime))'];
                        [ProdGross1, ProdNet1,ProdHTL1,eHTL,Bpico1,Bnano1,Bmicro1] = ...
                            getFunctions(u, sim.L(i,j,k,iTime));
                        conv = squeeze(dz(i,j,k));
                        ProdGross = ProdGross + ProdGross1*conv;
                        ProdNet = ProdNet + ProdNet1*conv;
                        ProdHTL = ProdHTL +ProdHTL1*conv;
                        %eHTL = eHTL + eHTL1/length(sim.z);
                        Bpico = Bpico + Bpico*dz(i,j,k); % gC/m2
                        Bnano = Bnano + Bnano1*dz(i,j,k);
                        Bmicro = Bmicro + Bmicro1*dz(i,j,k);
                    end
                end
                sim.ProdGross(i,j,iTime) = ProdGross;
                sim.ProdNet(i,j,iTime) = ProdNet;
                sim.ProdHTL(i,j,iTime) = ProdHTL;
                %sim.eHTL(i,j,iTime) = eHTL;
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
    function tot = calcTotal(u)
        tot = sum(u(ix(:)).*dv(ix(:))); % mug/day
    end

for i = 1:length(sim.t)
    sim.Ntotal(i) = calcTotal(sim.N(:,:,:,i));
    sim.DOCtotal(i) = calcTotal(sim.DOC(:,:,:,i)); % mugC
    sim.Btotal(i) = 0;
    for j = 1:sim.p.nGrid
        sim.Btotal(i) = sim.Btotal(i) + calcTotal(sim.B(:,:,:,j,i)); % mugC
    end
    sim.ProdNetTotal(i) = calcTotal(sim.ProdNet(:,:,i)); % mugC/day
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
% Annual global totals:
%
if (sim.t(end)>=365)
    sim.ProdNetAnnual = zeros(length(sim.x), length(sim.y), floor(sim.t(end)/365));
    for i = 1:sim.t(end)/365
        ixTime = sim.t>365*(i-1) & sim.t<=365*i;
        sim.ProdNetTotalAnnual(i) = mean(sim.ProdNetTotal(ixTime))/length(ixTime); % mugC/l/day
        sim.ProdNetAnnual(:,:,i) = mean(sim.ProdNet(:,:,ixTime),3)/length(ixTime); % gC/m2/yr
        sim.ProdHTLAnnual(:,:,i) = mean(sim.ProdHTL(:,:,ixTime),3)/length(ixTime); % gC/m2/yr
    end
    sim.ProdNetTotalAnnual = sim.ProdNetTotalAnnual/1000/1000/1000/1000; % GTon carbon/year
end

end