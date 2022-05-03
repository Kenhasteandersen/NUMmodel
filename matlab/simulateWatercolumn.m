%
% Simulate a single water column from the global transport matrix.
% Conservation is enforced rather crudely.
%
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% (choose MITgcm_ECCO_v4), and be put into the location 'NUMmodel/TMs'
%
% Input:
%  p: parameter structure from parametersGlobal
%  lat, lon: latitude and longitude
%  sim: (optional) simulation to use for initial conditions
%
% Input options:
%  bCalcAnnualAverages: increases the simulation time by a factor 2-3
%  bExtractcolumn: (logical) forces the extraction of a water column from
%            the transport matrix. Only used when the core code is changed.
%  dayFixed: if non-zero then the run is done with fixed conditions at the
%            specified day.
%
% Output:
%  sim: structure with simulation results
%
function sim = simulateWatercolumn(p, lat, lon, sim, options)

arguments
    p struct;
    lat = 0;
    lon = 0;
    sim struct = [];
    options.bExtractcolumn logical = false; % Extract the watercolumn even though a saved one exists
    options.bRecalcLight logical = false; % Recalc the light (different from the extracted watercolumn)
    options.dayFixed double = 0;
end
%
% Get the watercolumn parameters if they are not already set:
%
if ~isfield(p,'nameModel')
    p = parametersWatercolumn(p);
end

ixN = p.idxN;
ixDOC = p.idxDOC;

bSilicate = false;
if isfield(p,'idxSi')
    ixSi = p.idxSi;
    bSilicate = true;
end
ixB = p.idxB:p.n;
 
%Tbc = [];

disp('Preparing simulation')
%
% Set path:
%
path = fileparts(mfilename('fullpath'));
addpath(strcat(path,'/Transport matrix'));

% ---------------------------------------
% Find the indices that corresponds to the selected water column:
% ---------------------------------------
if isempty(sim)
    sim = load(p.pathGrid,'x','y','z','dznom','bathy'); % Load grid
    %sim.x = lon;
    %sim.y = lat;
end
load(p.pathBoxes, 'nb', 'Ybox', 'Zbox');
idx = calcGlobalWatercolumn(lat, lon, sim); % Find the indices into matrix
xx = matrixToGrid((1:nb)', [], p.pathBoxes, p.pathGrid); % Find the indices into the grid
idxGrid = squeeze(xx(idx.x, idx.y, idx.z));
nGrid = length(idxGrid);

if nGrid==0
    error('Selected latitude and longitude is on land.')
end
%
% Check if a file with the water column exists; if not extract from TM and
% save:
%
sFile = sprintf('Watercolumn_%s_lat%03i_lon%02i.mat',p.TMname,lat,lon);
if exist(sFile,'file') && ~options.bExtractcolumn
    load(sFile);
end

if ~exist('versionTMcolumn')
    versionTMcolumn=0;
end
versionTMcolumnCurrent = 2; % Current version of the water column
if (versionTMcolumn~=versionTMcolumnCurrent) || options.bExtractcolumn  % Extract water column if the loaded one is too old
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
        %disp(strcat(p.pathMatrix, sprintf('%02i.mat',month+1)));

        %AexpM(month,:,:) = full(function_convert_TM_positive(Aexp(idxGrid,idxGrid)));
        AimpM(month,:,:) = full(function_convert_TM_positive(matrix.Aimp(idxGrid,idxGrid)));

        % Preparing for timestepping. 43200s.
        temp = load(p.pathGrid,'deltaT');
        %AexpM(month,:,:) = Ix(idxGrid,idxGrid) + squeeze((12*60*60)*AexpM(month,:,:));
        AimpM(month,:,:) = squeeze(AimpM(month,:,:))^(p.dtTransport*24*60*60/temp.deltaT);
        fprintf('.');
    end
    %
    % Load temperature:
    %
    load(p.pathTemp, 'Tbc');
    Tmat = zeros(nb,12);
    for i = 1:12
        Tmat(:,i) = gridToMatrix(Tbc(:,:,:,i), [], p.pathBoxes, p.pathGrid);
    end
    Tmat = Tmat(idxGrid,:); % Use only the specific water column
    %
    % Calc Light:
    %
    L0 = zeros(nGrid,730);
    for i = 1:730
        zup = sim.z - 0.5*sim.dznom; % Top of a box
        zup = zup(1:length(idxGrid));
        dz = sim.dznom(1:length(idxGrid));
        Lup = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid),i/2,1).*exp(-p.kw*zup);
        L0(:,i) = Lup.*(1-exp(-p.kw*dz))./(p.kw*dz);
    end

    fprintf('\n');
    save(sFile,'AimpM','Tmat','L0','versionTMcolumn');
end

if options.bRecalcLight
    %
    % Calc Light:
    %
    L0 = zeros(nGrid,730);
    for i = 1:730
        zup = sim.z - 0.5*sim.dznom; % Top of a box
        zup = zup(1:length(idxGrid));
        dz = sim.dznom(1:length(idxGrid));
        Lup = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid),i/2,1).*exp(-p.kw*zup);
        L0(:,i) = Lup.*(1-exp(-p.kw*dz))./(p.kw*dz);
        %Lold = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid),i/2,1).*exp(-p.kw*Zbox(idxGrid));
    end
end
% Get sinking matrix:
[Asink,p] = calcSinkingMatrix(p, sim, nGrid);
% ---------------------------------------
% Initialize run:
% ---------------------------------------
simtime = p.tEnd/p.dtTransport; %simulation time in TM time steps

% Preparing timestepping
mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
%
% Initial conditions:
%
if isfield(sim,'B')
    disp('Starting from previous simulation.');
    u(:,ixN) = sim.N(:,end);
    u(:, ixDOC) = sim.DOC(:,end);
    if bSilicate
        u(:, ixSi) = sim.Si(:,end);
    end
    for i = 1:p.n -p.idxB+1
        u(:, ixB(i)) = sim.B(:,i,end);
    end
else
    if exist(strcat(p.pathN0,'.mat'),'file')
        load(p.pathN0, 'N');
        u(:, ixN) = gridToMatrix(N, [], p.pathBoxes, p.pathGrid);
    else
        u(:, ixN) = p.u0(p.idxN)*ones(nb,1);
    end
    u(:, ixDOC) = zeros(nb,1) + p.u0(ixDOC);
    if bSilicate
        u(:, ixSi) = zeros(nb,1) + p.u0(ixSi);
    end
    u(:, ixB) = ones(nb,1)*p.u0(ixB);
    u = u(idxGrid,:); % Use only the specific water column
    p.u0(ixN) = u(nGrid,ixN); % Use the nitrogen concentration in the last grid cell as BC
end
%
% Matrices for saving the solution:
%
iSave = 0;
nSave = floor(p.tEnd/p.tSave) + sign(mod(p.tEnd,p.tSave));

sim.N = zeros(length(idx.z),nSave);
if bSilicate
    sim.Si = sim.N;
end
sim.DOC = sim.N;
sim.B = zeros(length(idx.z), p.n-p.idxB+1, nSave);
sim.L = sim.N;
sim.T = sim.N;
sim.Nloss = zeros(1,nSave);
sim.NlossHTL = sim.Nloss;
tSave = [];
%
% Matrices for annual averages:
%
%if bCalcAnnualAverages
%   sim.ProdGrossAnnual = zeros( nb,1 );
%   sim.ProdNetAnnual = zeros( nb,1 );
%   sim.ProdHTLAnnual = zeros( nb,1 );
%   sim.BpicoAnnualMean = zeros( nb,1 );
%   sim.BnanoAnnualMean = zeros( nb,1 );
%   sim.BmicroAnnualMean = zeros( nb,1 );
%end

% ---------------------------------------
% Run transport matrix simulation
% ---------------------------------------
disp('Starting simulation')
tic
for i = 1:simtime
    %
    % Find the iteration time
    %
    if (options.dayFixed ~= 0)
        iTime = options.dayFixed / p.dtTransport;
    else
        iTime = i;
    end
    %
    % Test for time to change monthly temperature
    %
    if (ismember(mod(iTime,365/p.dtTransport), 1+2*cumsum(mon)) || i==1)
        %fprintf('t = %u days\n',floor(i/2))
        % Set monthly mean temperature
        month = find(1+2*cumsum(mon) >= mod(iTime,365/p.dtTransport),1);
        T = Tmat(:,month);
        %month = mod(month + 1, 12);
    end
    %
    % Run Euler time step for half a day:
    %
    L = L0(:,mod(iTime,365*2)+1);
    dt = p.dt;
    dtTransport = p.dtTransport;
    n = p.n;
    %if ~isempty(gcp('nocreate'))
    %    parfor k = 1:nGrid
    %        u(k,:) = calllib(loadNUMmodelLibrary(), 'f_simulateeuler', ...
    %            int32(n), u(k,:), L(k), T(k), 0.5, dt);
    %    end
    %else
    for k = 1:nGrid
        u(k,:) = calllib(loadNUMmodelLibrary(), 'f_simulateeuler', ...
            u(k,:),L(k), T(k), dtTransport, dt);
    end
    %end

    if any(isnan(u))
        warning('NaNs after running current time step');
        keyboard
    end
    %
    % Transport
    %
    u =  squeeze(AimpM(month,:,:)) * u; % Vertical diffusion
    % Sinking:
    for j = p.idxSinking
        u(:,j) = squeeze(Asink(j,:,:)) * u(:,j);
    end
    % Bottom BC for nutrients:
    u(end, p.idxN) = u(end, p.idxN) +  p.dtTransport* ...
        p.DiffBottom/sim.dznom(nGrid)*(p.u0(p.idxN)-u(end,p.idxN));
    %
    % Enforce minimum concentration
    %
    for k = 1:nGrid
        u(k,u(k,:)<p.umin) = p.umin(u(k,:)<p.umin);
    end
    %
    % Save timeseries in grid format
    %
    if ((mod(i/2,p.tSave) < mod((i-1)/2,p.tSave)) || (i==simtime))
        iSave = iSave + 1;
        sim.N(:,iSave) = u(:,ixN);
        sim.DOC(:,iSave) = u(:,ixDOC);
        if bSilicate
            sim.Si(:,iSave) = u(:,ixSi);
        end
        for j = 1:p.n-p.idxB+1
            sim.B(:,j,iSave) = u(:,ixB(j));
        end
        sim.L(:,iSave) = L;
        sim.T(:,iSave) = T;
        % Loss to HTL and remin2:
        for j = 1:length(nGrid)
            rates = getRates(p,u(j,:),L(j),T(j));
            % Note: half of the HTL loss is routed directly back to N if we
            % don't have POM:
            remin2 = 0.5;
            reminHTL = 0.5;
            rhoCN=5.68;
            if ~sum(ismember(p.typeGroups,100))
                sim.NlossHTL(iSave) = sim.NlossHTL(iSave) + ...
                    reminHTL*sum(rates.mortHTL.*u(j,p.idxB:end)')/1000*sim.dznom(j)/rhoCN; % % HTL losses:gN/m2/day
                sim.Nloss(iSave) = sim.Nloss(iSave) + remin2*sum(rates.mort2.*u(j,p.idxB:end)')/1000*sim.dznom(j)/rhoCN; % remin2 losses
            end
        end

        tSave = [tSave, i*p.dtTransport];
    end
    %
    % Update annual averages:
    %
    %     if bCalcAnnualAverages
    %         for k = 1:nb
    %             [ProdGross1, ProdNet1,ProdHTL1,eHTL,Bpico1,Bnano1,Bmicro1] = ...
    %                             getFunctions(u(k,:), L(k));
    %             sim.ProdGrossAnnual(k) = sim.ProdGrossAnnual(k) + ProdGross1/(p.tEnd*2);
    %             sim.ProdNetAnnual(k) = sim.ProdNetAnnual(k) + ProdNet1/(p.tEnd*2);
    %             sim.ProdHTLAnnual(k) = sim.ProdHTLAnnual(k) + ProdHTL1/(p.tEnd*2);
    %             sim.BpicoAnnualMean(k) = sim.BpicoAnnualMean(k) + Bpico1/(p.tEnd*2*365);
    %             sim.BnanoAnnualMean(k) = sim.BnanoAnnualMean(k) + Bnano1/(p.tEnd*2*365);
    %             sim.BmicroAnnualMean(k) = sim.BmicroAnnualMean(k) + Bmicro1/(p.tEnd*2*365);
    %         end
    %     end

end
time = toc;
fprintf('Solving time: %2u:%02u:%02u\n', ...
    [floor(time/3600), mod(floor(time/60),60), floor(mod(time,60))]);
% ---------------------------------------
% Put results into sim structure:
% ---------------------------------------
sim.t = tSave; % days where solution was saved
sim.p = p;
%sim.Ntot = calcGlobalN(sim);
sim.B(sim.B<0) = 0.;
sim.DOC(sim.DOC<0) = 0.;
sim.z = sim.z(1:length(idx.z));
sim.dznom = sim.dznom(1:length(idx.z));
sim.lat = lat;
sim.lon = lon;

sim.Ntot = (sum(sim.N.*(sim.dznom*ones(1,length(sim.t)))) + ... % gN/m2 in dissolved phase
    sum(squeeze(sum(sim.B,2)).*(sim.dznom*ones(1,length(sim.t))))/5.68)/1000; % gN/m2 in biomass
sim.Nprod = p.DiffBottom*(p.u0(p.idxN)-sim.N(end,:))/1000; % Diffusion in from the bottom; gN/m2/day
% if bCalcAnnualAverages
%     tmp = single(matrixToGrid(sim.ProdGrossAnnual, [], p.pathBoxes, p.pathGrid));
%     sim.ProdGrossAnnual = squeeze(tmp(:,:,1));
%     tmp = single(matrixToGrid(sim.ProdNetAnnual, [], p.pathBoxes, p.pathGrid));
%     sim.ProdNetAnnual = squeeze(tmp(:,:,1));
%     tmp = single(matrixToGrid(sim.ProdHTLAnnual, [], p.pathBoxes, p.pathGrid));
%     sim.ProdHTLAnnual = squeeze(tmp(:,:,1));
%     tmp = single(matrixToGrid(sim.BpicoAnnualMean, [], p.pathBoxes, p.pathGrid));
%     sim.BpicoAnnualMean = squeeze(tmp(:,:,1));
%     tmp = single(matrixToGrid(sim.BnanoAnnualMean, [], p.pathBoxes, p.pathGrid));
%     sim.BnanoAnnualMean = squeeze(tmp(:,:,1));
%     tmp = single(matrixToGrid(sim.BmicroAnnualMean, [], p.pathBoxes, p.pathGrid));
%     sim.BmicroAnnualMean = squeeze(tmp(:,:,1));
% end
end
