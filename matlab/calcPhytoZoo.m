%
% Return the biomass of a each group that can be associated with
% phototrophy ("phytoplankton"), osmotrophy ("bacteria"), and phagotrophy
% ("zooplankton").
%
% In:
%  p - parameters
%  u - state vector
%  L - Light
%  T - temperature
%
% Out:
%  Bphyto, Bzoo, Bbacteria - Biomass in mugC/l
%
function [Bphyto, Bzoo, Bbacteria] = calcPhytoZoo(p, u, L, T, sLibname)
arguments
    p struct;
    u double;
    L double;
    T double;
    sLibname = loadNUMmodelLibrary();
end

rates = getRates(p, u, L, T, sLibname);
rates.jLreal( isnan(rates.jLreal) ) = 0;
jCarbon = rates.jDOC + rates.jFreal + rates.jLreal + 1e-100;

B = u(p.idxB:end)';
for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup)) - p.idxB + 1;

    if p.typeGroups(iGroup)<10 % For unicellular groups
        Bphyto(iGroup) = sum( B(ix).*rates.jLreal(ix)./jCarbon(ix) );
        Bbacteria(iGroup) = sum( B(ix).*rates.jDOC(ix)./jCarbon(ix) );
        Bzoo(iGroup) = sum( B(ix).*rates.jFreal(ix)./jCarbon(ix) );
    else
        Bphyto(iGroup) = 0;
        Bbacteria(iGroup) = 0;
        Bzoo(iGroup) = sum( B(ix) );
    end
end


