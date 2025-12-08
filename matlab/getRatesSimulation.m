%
% Get the rates from a simulation at a given time, depth and position
%
% In:
%   sim - simulation structure
%   iTime - Time step
%   iDepth - The depth stratum (for watercolumn and global simulation)
%   iLon, iLat - Longitude and latitude indices (for global simulation)
%
% Out:
%   rates - The rates for each size groups (1/day)
%   B - The biomasses
%
function [rates, B] = getRatesSimulation(sim, iTime, iDepth, iLon, iLat)

sLibName = loadNUMmodelLibrary();

switch sim.p.nameModel

    case 'chemostat'
        if length(sim.L)>1
            L = sim.L(iTime);
        else
            L = sim.L;
        end
        if length(sim.T)>1
            T = sim.T(iTime);
        else
            T = sim.T;
        end

        B = sim.B(iTime,:);
        if sim.p.nNutrients==3
            u = [sim.N(iTime), sim.DOC(iTime),sim.Si(iTime), B];
        else
            u = [sim.N(iTime), sim.DOC(iTime), B];
        end
        rates = getRates(sim.p, u, L, T, sLibName);
        B = B';

    case 'watercolumn'
        % Get the functions per volume at each depth and time:
        B = squeeze(sim.B(iTime,iDepth,:));
        if sim.p.nNutrients==3
            u = [squeeze(sim.N(iTime,iDepth)), ...
                squeeze(sim.DOC(iTime,iDepth)), ...
                squeeze(sim.Si(iTime,iDepth)), ...
                B'];
        else
            u = [squeeze(sim.N(iTime,iDepth)), ...
                squeeze(sim.DOC(iTime,iDepth)), ...
                B'];
        end
        rates = getRates(sim.p, u, sim.L(iTime,iDepth), sim.T(iTime,iDepth), ...
            sLibName);

    case 'global'
        B = squeeze(sim.B(iTime,iLon, iLat, iDepth,:));
        if isfield(sim.p,'idxSi')
            u = [squeeze(sim.N(iTime,iLon, iLat, iDepth)), ...
                squeeze(sim.DOC(iTime,iLon, iLat, iDepth)), ...
                squeeze(sim.Si(iTime,iLon, iLat, iDepth)), ...
                B'];
        else
            u = [squeeze(sim.N(iTime,iLon, iLat, iDepth)), ...
                squeeze(sim.DOC(iTime,iLon, iLat, iDepth)), ...
                B'];
        end
        rates = getRates(sim.p, u, sim.L(iTime,iLon, iLat, iDepth), ...
            sim.T(iTime,iLon, iLat, iDepth), ...
            sLibName);

    otherwise
        disp('Model unknown; rates not calculated');

end

