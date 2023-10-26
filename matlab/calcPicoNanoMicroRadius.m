%
% Returns the biomass of pico, nano and micro plankton in mugC/liter.
% Only works with generalists.
% Example: call calcPicoNanoMicroRadius(sim.B(end,:), sim.p.pGeneralists, groups);
% or B=squeeze(mean(simNUM.B(:,1,:),1));

% options: select group, all groups,...
function Bpnm = calcPicoNanoMicroRadius(B,p,groups)
% arguments
%     p struct;
%     B = [];
%     groups = false;
% end

% Calc radius
r = calcRadiusGroups(p);
% Masses for pico, nano, and micro plankton
r0 = 0;
r2 = 1;
r20 = 10;
r200 = 100;

Bpico = calcBiomassRangeRadius(B, r, r0,r2);
Bnano = calcBiomassRangeRadius(B, r, r2,r20);
Bmicro = calcBiomassRangeRadius(B, r, r20,r200);

if groups==false
Bpnm = [Bpico, Bnano, Bmicro];
end

if groups==true
    BpicoGroup=zeros(1,p.nGroups);
    BnanoGroup=zeros(1,p.nGroups);
    BmicroGroup=zeros(1,p.nGroups);
    Bpnm = zeros(p.nGroups,3);
   for iGroup = 1:p.nGroups
         ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
       % check if there is not POM
       if p.typeGroups(iGroup)~=100
        BpicoGroup(iGroup) = calcBiomassRangeRadius(B(ix), r(ix), r0,r2);
        BnanoGroup(iGroup) = calcBiomassRangeRadius(B(ix), r(ix), r2,r20);
        BmicroGroup(iGroup) = calcBiomassRangeRadius(B(ix), r(ix), r20,r200);
       else
           if ( r(ix)>=r0 && r(ix)<r2)
               BpicoGroup(iGroup)= B(ix);
           end
           if ( r(ix)>=r2 && r(ix)<r20)
               BnanoGroup(iGroup)= B(ix);
           end
           if ( r(ix)>=r20 && r(ix)<r200)
               BmicroGroup(iGroup)= B(ix);
           end
       end
        Bpnm(iGroup,:)=[BpicoGroup(iGroup), BnanoGroup(iGroup), BmicroGroup(iGroup)];
    end
 end
end

% plotGlobalBiomass(sim, sProjection)