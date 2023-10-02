%
% Returns the biomass of pico, nano and micro plankton in mugC/liter.
% Bph = calcPhyto(sim,simR)
% example of Bnum: Bnum=squeeze(mean(simNUM.B(:,1,:),1));
%
% options: select group, all groups,...
function Bpnm = calcPhytoPicoNanoMicroRadius(Bph,p,sim)
% arguments
%     B = [];
%     p struct;
%     sim struct;
% end


[~,idxU]=find(sim.p.typeGroups<10);
ixUPlankton = sim.p.ixStart(idxU(1)):sim.p.ixEnd(idxU(end));
%%


% Calc radius
r_all = calcRadiusGroups(p);
r = r_all(ixUPlankton); 
% Masses for pico, nano, and micro plankton
r0 = 0;
r2 = 1;
r20 = 10;
r200 = 100;

Bpico = calcBiomassRangeRadius(Bph, r, r0,r2);
Bnano = calcBiomassRangeRadius(Bph, r, r2,r20);
Bmicro = calcBiomassRangeRadius(Bph, r, r20,r200);

Bpnm = [Bpico, Bnano, Bmicro];