%
% Returns the biomass of organisms that get photosynthesize more than 60%
% simR=getSimRates(sim);
function Bph = calcPhyto(sim,simR)

% arguments
%     p struct;
% 
% end


sLibName = loadNUMmodelLibrary();
ixTime = find(sim.t>(max(sim.t)-365)); %nTime = length(sim.t(sim.t >= max(sim.t-365))); % Just do the last year
nZ = length(sim.z);
[~,idxU]=find(sim.p.typeGroups<10);
ixUPlankton = sim.p.ixStart(idxU(1)):sim.p.ixEnd(idxU(end));
Bph=sim.B(:,:,ixUPlankton);

  for iTime = ixTime
      for k = 1:nZ
         Bphyto_frac(iTime,k,:)=calcBphytoFrac( simR.jLreal(iTime,k,:),...
         simR.jF(iTime,k,:),simR.jDOC(iTime,k,:),ixUPlankton);
      end
  end

for iTime = ixTime
            for k = 1:nZ
                for i=1:length(ixUPlankton)
                    if Bphyto_frac(iTime,k,i)<0.6
                        Bph(iTime,k,i)=0;
                    end
                end

            end
end



function Bphyto_frac = calcBphytoFrac(jLreal,jF,jDOC,ix) %check units!
        Bphyto_frac = jLreal(ix)./(jLreal(ix)+jF(ix)+jDOC(ix));
end

end