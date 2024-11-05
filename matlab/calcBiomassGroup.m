%
% Calculate the total biomass in group per area (generalists, active copepods, etc.)
%
% In:
%  sim - simulation output structure
%  typeGroup - The type of the group (see parametersAddGroup for definition
%              of groups) 
%  bAverageTime - whether to average over the last year (default true)
% 
% Out:
%  B - biomass in mugC/area
%
function B = calcBiomassGroup(sim, typeGroup, bAverageTime)

arguments
    sim struct;
    typeGroup int32;
    bAverageTime = true;
end
%
% Sum over group:
% 
ixGroup = find( sim.p.typeGroups==typeGroup );
if isempty(ixGroup)
    error('There is no %i group in the simulation. Check sim.p.typeGroups\n',typeGroup);
end
ix = [];
for i = 1:length(ixGroup)
    ix = [ix (sim.p.ixStart(ixGroup(i)):sim.p.ixEnd(ixGroup(i)))-sim.p.idxB+1];
end
%
% Prepare time averaging:
%
if bAverageTime
    ixTime = find(sim.t > sim.t(end)-365);
else
    ixTime = length(sim.t);
end
%
% Integrate over depth:
%
switch sim.p.nameModel
    
    case 'chemostat'
        B = sum( sim.B(:, ix),2 ) * sim.p.widthProductiveLayer;
        B = mean( B(ixTime,:),1 );

    case 'watercolumn'
        B = squeeze( sum( sim.B .* reshape(sim.dznom,1,numel(sim.dznom),1),2) );
        B = sum( B(:, ix),2 );
        B = mean( B(ixTime,:),1 );

    case 'global'
        B = sum( sim.B(:,:,:,:,ix),5 ); % Sum the group
        B = squeeze( calcIntegrateGlobal(sim,B,bAverageTime) );
end





