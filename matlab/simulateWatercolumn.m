%
% Simulate a single water column from the global transport matrix.
% Conservation is enforced rather crudely.
%
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% (choose MITgcm_ECCO), and put it into the location 'NUMmodel/TMs/'
%
% Input:
%  p: parameter structure
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
    %options.bRecalcLight logical = false; % Recalc the light (different from the extracted watercolumn)
    options.dayFixed double = 0;
end
disp('Preparing simulation')

%
% Get the watercolumn parameters if they are not already set:
%
if ~isfield(p,'nameModel')
    p = parametersWatercolumn(p);
end
%
% Define some shorthands:
%
ixN = p.idxN;
ixDOC = p.idxDOC;

bSilicate = false;
if isfield(p,'idxSi')
    ixSi = p.idxSi;
    bSilicate = true;
end
ixB = p.idxB:p.n;

S = inputRead;
reminHTL = S.input_general.fracHTL_to_N;
rhoCN = S.input_general.rhoCN;
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
    fprintf('-> Extracting water column from transport matrix\n')
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
    %L0 = calclight;

    fprintf('\n');
    save(sFile,'AimpM','Tmat','versionTMcolumn');
end

L0 = calclight;

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
    u(:,ixN) = sim.N(end,:);
    u(:, ixDOC) = sim.DOC(end,:);
    if bSilicate
        u(:, ixSi) = sim.Si(end,:);
    end
    for i = 1:p.n -p.idxB+1
        u(:, ixB(i)) = sim.B(end,:,i);
    end
else
    if exist(strcat(p.pathN0,'.mat'),'file')
        load(p.pathN0, 'N');
        u(:, ixN) = gridToMatrix(N, [], p.pathBoxes, p.pathGrid);
    else
        u(:, ixN) = p.u0(p.idxN)*ones(nb,1);
    end

    if bSilicate
        if exist(strcat(p.pathSi0,'.mat'),'file')
            load(p.pathSi0, 'Si');
            u(:, ixSi) = gridToMatrix(Si, [], p.pathBoxes, p.pathGrid);
        else
            u(:, ixSi) = zeros(nb,1) + p.u0(ixSi);
        end
    end

    u(:, ixDOC) = zeros(nb,1) + p.u0(ixDOC);
    u(:, ixB) = ones(nb,1)*p.u0(ixB);
    u = u(idxGrid,:); % Use only the specific water column
end
%p.u0(ixN) = u(nGrid,ixN); % Use the nitrogen concentration in the last grid cell as BC
if bSilicate
    p.u0(ixSi) = u(nGrid,ixSi);
end
%
% Set BCvalue:
%
BCvalue = p.BCvalue;
%if size(BCvalue,1)==1
%    BCvalue = ones(length(ixBottom),1)*BCvalue;
%end
% If BCvalue == -1 then use the bottom value from the initial conditions:
for i = 1:length(BCvalue)
    if BCvalue(i)==-1
        BCvalue(i) = u(end,i)';
    end
end
%
% Matrices for saving the solution:
%
iSave = 0;
nSave = floor(p.tEnd/p.tSave) + sign(mod(p.tEnd,p.tSave));

sim.N = zeros(nSave, length(idx.z));
if bSilicate
    sim.Si = sim.N;
end
sim.DOC = sim.N;
sim.B = zeros(nSave, length(idx.z), p.n-p.idxB+1);
sim.L = sim.N;
sim.T = sim.N;
sim.Nloss = zeros(nSave,1);
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
sLibName = loadNUMmodelLibrary();
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
    L = L0(:,mod(iTime,365/p.dtTransport)+1);
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
        u(k,:) = calllib(sLibName, 'f_simulateeuler', ...
            u(k,:),L(k), T(k), dtTransport, dt);
    end

    %if sum(u(:)<0)
    %    fprintf("Time step %i, #i\n", [i, sum(u(:)<0)]);
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
    u(end, 1:p.nNutrients) = u(end, 1:p.nNutrients) +  ...
        p.BCmixing(1:p.nNutrients)*p.dtTransport .* ...
        ( BCvalue(1:p.nNutrients) - u(end,1:p.nNutrients) );

    %u(end, p.idxN) = u(end, p.idxN) +  p.dtTransport* ...
    %    p.DiffBottom/sim.dznom(nGrid)*(p.u0(p.idxN)-u(end,p.idxN));

    %if bSilicate
    %    u(end, p.idxSi) = u(end, p.idxSi) +  p.dtTransport* ...
    %    p.DiffBottom/sim.dznom(nGrid)*(p.u0(p.idxSi)-u(end,p.idxSi));
    %end

    %
    % Enforce minimum concentration
    %
    %for k = 1:nGrid
    %    u(k,u(k,:)<p.umin) = p.umin(u(k,:)<p.umin);
    %end
    %
    % Save timeseries in grid format
    %
    if ((floor(i*(p.dtTransport/p.tSave)) > floor((i-1)*(p.dtTransport/p.tSave))) || (i==simtime))
        iSave = iSave + 1;
        sim.N(iSave,:) = u(:,ixN);
        sim.DOC(iSave,:) = u(:,ixDOC);
        if bSilicate
            sim.Si(iSave,:) = u(:,ixSi);
        end
        for j = 1:p.n-p.idxB+1
            sim.B(iSave,:,j) = u(:,ixB(j));
        end
        sim.L(iSave,:) = L;
        sim.T(iSave,:) = T;
        % Loss to HTL and POM:
        for j = 1:nGrid
            rates = getRates(p,u(j,:),L(j),T(j));
            % Note: half of the HTL loss is routed directly back to N if we
            % don't have POM:
            if ~sum(ismember(p.typeGroups,100))
                sim.NlossHTL(iSave) = sim.NlossHTL(iSave) + ...
                    (1-reminHTL)*sum(rates.mortHTL.*u(j,p.idxB:end)')/1000*sim.dznom(j)/rhoCN; % % HTL losses:gN/m2/day
                sim.Nloss(iSave) = sim.Nloss(iSave) + ...
                    sum(rates.jPOM.*u(j,p.idxB:end)')/1000*sim.dznom(j)/rhoCN;
                %remin2*sum(rates.mort2.*u(j,p.idxB:end)')/1000*sim.dznom(j)/rhoCN; % remin2 losses
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

sim.Ntot = (sum(sim.N'.*(sim.dznom*ones(1,length(sim.t)))) + ... % gN/m2 in dissolved phase
    sum(squeeze(sum(sim.B,3))'.*(sim.dznom*ones(1,length(sim.t))))/rhoCN)/1000; % gN/m2 in biomass
sim.Nprod = p.BCmixing(p.idxN)*(p.BCvalue(p.idxN)-sim.N(:,end))*sim.dznom(end)/1000; % Diffusion in from the bottom; gN/m2/day
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



%
% Setup matrix for sinking. Uses an implicit first-order upwind scheme
%
    function [Asink,p] = calcSinkingMatrix(p, sim, nGrid)
        %
        % Get sinking velocities from libNUMmodel:
        %
        p.velocity = 0*p.m;
        p.velocity = calllib(loadNUMmodelLibrary(), 'f_getsinking', p.velocity);
        p.idxSinking = find(p.velocity ~= 0); % Find indices of groups with sinking
        %
        % Set up the matrix:
        %
        Asink = zeros(p.n,nGrid,nGrid);
        for i = 1:nGrid
            k = p.velocity*p.dtTransport/sim.dznom(i);
            Asink(:,i,i) = 1+k;
            if (i ~= 1)
                Asink(:,i,i-1) = -k;
            end
        end
        % Bottom BC:
        if p.BC_POMclosed
            Asink(:,end,end) = 1;
        end
        %
        % Invert matrix to make it ready for use:
        %
        for i = 1:p.n
            Asink(i,:,:) = inv(squeeze(Asink(i,:,:)));
        end
    end

    function L0 = calclight()
        L0 = zeros(nGrid,365/p.dtTransport);

        if p.bUse_parday_light
            if exist(p.pathPARday,'file')
                load(p.pathPARday,'parday');
            else
                error('PARday file does not exist. Set p.bUse_parday_light = false');
            end
        end

        for i = 1:365/p.dtTransport
            % zup = sim.z - 0.5*sim.dznom; % Top of a box
            % zup = zup(1:length(idxGrid));
            % dz = sim.dznom(1:length(idxGrid));
            % if p.bUse_parday_light
            %     Lup = 1e6*parday(idx.x,idx.y,1,i)/(24*60*60).*exp(-p.kw*zup);
            % else
            %     Lup = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid),i/2,1).*exp(-p.kw*zup);
            % end
            % L0(:,i) = Lup.*(1-exp(-p.kw*dz))./(p.kw*dz);
            if p.bUse_parday_light
                Lup = 1e6*parday(idx.x,idx.y,1,i)/(24*60*60);
            else
                Lup = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid(1)),i*p.dtTransport,1);
            end
            L0(:,i) = Lup*exp(-p.kw*sim.z(1:nGrid));
        end

    end

end

