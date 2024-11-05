%
% Set the parameters for a chemostat simulation.
%
% In:
%  p - parameter structure from a call to a setupXX function
%  seasonalOptions - parameter structure that choose if the chemostat will
%                    have or not seasonalities. Several types are possible:
%                    (author: Cécile Decker)
%  seasonalOptions.constantValues - set to 0.5 and 100 by default (mixing
%                                   rate and light are constant)
%  seasonalOptions.lat_lon - set to NaN by default (mixing rate and light 
%                            are constant, respectively 0.5 and 100). For a 
%                            mixing rate and light that correspond to a 
%                            particular location and are calculated 
%                            seasonaly enter an array of two integer (lat 
%                            and lon).
%  seasonalOptions.seasonalAmplitude - set to 0 by default (mixing rate and
%                                      light are constant, respectively 0.5
%                                      and 100). You can choose a double
%                                      between 0 and 1 the mixing rate and
%                                      light will then correspond to a 
%                                      weighted sum between the previous
%                                      constant and a simplified mean
%                                      version of a seasonal zone.
%
% Out:
%  p - parameter structure with chemostat fields added.
%
function p = parametersChemostat(p, seasonalOptions)

arguments
    p struct = parametersChemostat(setupGeneralistsOnly);
    seasonalOptions.constantValues = [0.5 100]; % lower the mixing rate
    seasonalOptions.lat_lon = NaN;
    seasonalOptions.seasonalAmplitude = 0;
end

p.nameModel = 'chemostat';

p.d = 0.1;  % Default mixing rate (1/days)

p.widthProductiveLayer = 20; % (meters) only used for calcFunctions
p.tEnd = 365;  % Time to run in days
p.tSave = 1;

p.uDeep(1:p.idxB-1) = p.u0(1:p.idxB-1); % Nutrients of the layer below the chemostat layer
p.u0(1) = 1; % Nitrogen concentration

%
% Set minimum concentrations:
%
%p.umin = 0*p.u0;
%for iGroup = 1:p.nGroups
%    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
%    if p.typeGroups(iGroup) < 10
%        p.umin(ix) = 1e-5*p.mDelta(ix(1))/p.m(ix(1)); % Minimum B concentration for unicellular groups
%    end
%    if p.typeGroups(iGroup)>=10 && p.typeGroups(iGroup)<100
%        p.umin(ix(1)) = 1e-5*p.mDelta(ix(1))/p.m(ix(1)); % Send in some nauplii in copepod groups
%    end
%end


%
% Light environment:
%
% Parameters used to calculate light if not using parday:
p.EinConv = 4.57; % conversion factor from W m^-2 to \mu mol s^-1 m^-2 (Thimijan & Heins 1983)
p.PARfrac = 0.4; % Fraction of light available as PAR. Source unknown
p.kw = 0.05; % Damping of light by water; m^-1
%
% Check that transport matrix files exist:
%
path = fileparts(mfilename('fullpath'));
addpath(strcat(path,'/Transport matrix'));
if ~exist(path,'file')
    error( 'Error: Cannot find transport matrix file: %s', path);
end

% Mixing rate and light calculation regarding the option choose
depthProductiveLayer = 50; % (meters) Depth of the productive layer. Used for seasonal calculations
if isnan(seasonalOptions.lat_lon)
    if seasonalOptions.seasonalAmplitude == 0
        p.d = 0.5;  % Mixing rate (1/days)
        p.L = 100; % constant default light
    else
        lambda = seasonalOptions.seasonalAmplitude;
        d = [];  % Mixing rate (1/days)
        L = [];  % Light
        %
        % In this option the mixing rate is a simplified seasonal version,
        % the idea was to choose the amplitude of this seasonality so we
        % are taking the barycenter between the previous constant value and
        % this new seasonal one using the option as the weight.
        %
        for t=1:365
            if t<125 | t>300
                d(t) = (1-lambda)*0.5+lambda*0.6;
            else
                d(t) = (1-lambda)*0.5+lambda*0.02;
            end
            if t<=60 | t>=280
                L(t) = (1-lambda)*100+lambda*5;
            else
                L(t) = (1-lambda)*100+lambda*35;
            end
        end
        p.d = d;
        p.L = L;
    end
else
    clear lat lon
    lat = seasonalOptions.lat_lon(1);
    lon = seasonalOptions.lat_lon(2);
    %
    % For this option the mixing rate is taken from the transport matrix of
    % this particular location (lat-lon)
    %
    p.TMname = 'MITgcm_ECCO';
    p.pathMatrix = strcat(path,'/../TMs/MITgcm_ECCO/Matrix1/TMs/matrix_nocorrection_');
    p.pathBoxes = strcat(path,'/../TMs/MITgcm_ECCO/Matrix1/Data/boxes.mat');
    p.pathGrid = strcat(path,'/../TMs/MITgcm_ECCO/grid.mat');
    p.dtTransport = 0.5; % The TM time step (in units of days)
    % Set path:
    path = fileparts(mfilename('fullpath'));
    addpath(strcat(path,'/Transport matrix'));
    grid = load(p.pathGrid,'x','y','z','dznom','bathy'); % Load grid
    load(p.pathBoxes, 'nb', 'Ybox', 'Zbox');
    idx = calcGlobalWatercolumn(lat, lon, grid); % Find the indices into matrix
    xx = matrixToGrid((1:nb)', [], p.pathBoxes, p.pathGrid); % Find the indices into the grid
    idxGrid = squeeze(xx(idx.x, idx.y, idx.z));
    nGrid = length(idxGrid);
    if nGrid==0
        error('Selected latitude and longitude is on land.')
    end
    sFile = sprintf('Watercolumn_%s_lat%03i_lon%02i.mat',p.TMname,lat,lon);
    if exist(sFile,'file')
        load(sFile);
    end
    if ~exist('versionTMcolumn')
        versionTMcolumn=0;
    end
    versionTMcolumnCurrent = 2; % Current version of the water column
    if (versionTMcolumn~=versionTMcolumnCurrent)  % Extract water column if the loaded one is too old
        versionTMcolumn = versionTMcolumnCurrent;
        fprintf('-> Extracting water column from transport matrix')
        %
        % Check that transport matrix files exist:
        %
        path = fileparts(mfilename('fullpath'));
        addpath(strcat(path,'/Transport matrix'));
        if ~exist(p.pathBoxes,'file')
            error( 'Error: Cannot find transport matrix file: %s', p.pathBoxes);
        end
        % Load TMs
        parfor month=1:12
            matrix = load(strcat(p.pathMatrix, sprintf('%02i.mat',month)),'Aimp');
            AimpM(month,:,:) = full(function_convert_TM_positive(matrix.Aimp(idxGrid,idxGrid)));
            temp = load(p.pathGrid,'deltaT');
            AimpM(month,:,:) = squeeze(AimpM(month,:,:))^(p.dtTransport*24*60*60/temp.deltaT);
            fprintf('.');
        end
    end
    z = grid.z;
    epsilon = 15;
    idx_chemo = find(z<(depthProductiveLayer+epsilon) & z>(depthProductiveLayer-epsilon));
    mon = [31 28 31 30 31 30 31 31 30 31 30 31];
    for iMonth = 1:12
        A = squeeze(AimpM(iMonth,:,:));
        j = sum(mon(1:iMonth-1));
        for k=1:mon(iMonth)
            d(k+j) = sum(A(idx_chemo,idx_chemo+1:end)); % sum all the layer below depth of the produtive layer
            %
            % Light function : we take 100m as the depth of the layer and a
            % layer of 20m width, because the chemostat is this all layer.
            %
            L(k+j) = (p.EinConv*p.PARfrac*daily_insolation(0,lat,k+j,1)*exp(-p.kw*depthProductiveLayer));%*(1-exp(-p.kw*p.widthProductiveLayer))/(p.kw*p.widthProductiveLayer);
        end
    end
    p.d = d;
    p.L = L;
end

p.seasonalOptions = seasonalOptions;

end
