%
% Global run using transport matrices
%
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% and be put into the location 'NUMmodel/TMs'
% Simulations currently works with:
%  - MITgcm_2.8deg (low resolution; runs on a laptop)
%  - MITgcm_ECCO (higher resolution; requires more memory).
%
% Input:
%  p: parameter structure from parametersGlobal
%  sim: (optional) simulation to use for initial conditions. Note the 
%            file with initial conditions for NUMmodel setup in the TM 
%            directory.
%  options.bCalcAnnualAverages: This option calculates ecosystem functions
%            (sim.ProdNet, sim.ProdNetAnnual etc) during the simulation, 
%            thus making them more accurate. It increases the simulation 
%            time of the last year.
%  options.bVerbose: Write out progress on the terminal
%  options.bNinit: Use only nutrient fields for initialization
%
% Output:
%  sim: structure with simulation results
%
function sim = simulateGlobal(p, sim, options)

arguments
    p struct;
    sim struct = [];
    options.bCalcAnnualAverages = false; % Whether to calculate annual averages
    options.bVerbose = true; % Whether to write output on the terminal
    options.bNinit = false;
end
%
% Get the global parameters if they are not already set:
%
if ~isfield(p,'nameModel')
    disp('Parameters not set. Using standard parameters for global run.');
    p = parametersGlobal(p);
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

if options.bVerbose
    disp('Preparing simulation')
end
%
% Check that files exist:
%
path = fileparts(mfilename('fullpath'));
addpath(strcat(path,'/Transport matrix'));

if ~exist(p.pathBoxes,'file')
    error( sprintf('Error: Cannot find transport matrix file: %s',...
        p.pathBoxes));
end
%
% Check that the time steps fit together. Everything below counts whole steps,
% and a step that does not divide evenly is caught nowhere else: it simply gives
% a run that is wrong. This is checked here rather than in parametersGlobal so
% that it also catches p.dtTransport or p.dt being overridden afterwards.
%
bWhole = @(x) abs(x-round(x)) < 1e-9;

if ~bWhole(p.dtTransport/p.dt)
    error(['p.dtTransport (%g days) must be a whole number of euler steps p.dt ' ...
           '(%g days), but it is %g of them. The library integrates ' ...
           'floor(dtTransport/dt) steps, so the biology would lag the transport.'], ...
           p.dtTransport, p.dt, p.dtTransport/p.dt);
end
if ~bWhole(365/p.dtTransport)
    error(['A year has to be a whole number of transport steps, but 365 days is ' ...
           '%g steps of p.dtTransport = %g days.'], 365/p.dtTransport, p.dtTransport);
end
%
% The implicit matrix is raised to the power dtTransport/deltaT below, and a
% sparse matrix can only be raised to a whole power:
%
gridDeltaT = load(p.pathGrid, 'deltaT');
if isfield(gridDeltaT,'deltaT')
    nPower = p.dtTransport*24*60*60/gridDeltaT.deltaT;
    if ~bWhole(nPower) || nPower < 1
        error(['p.dtTransport (%g days) is %g times the %g s time step of the ' ...
               '%s matrices. It has to be a whole multiple of it, because the ' ...
               'implicit matrix is raised to that power.'], ...
               p.dtTransport, nPower, gridDeltaT.deltaT, p.TMname);
    end
end
% ---------------------------------------
% Initialize run:
% ---------------------------------------
simtime = p.tEnd/p.dtTransport; % simulation time in transport steps
load(p.pathBoxes, 'nb', 'Ybox', 'Zbox');

% Preparing timestepping
Ix = speye(nb,nb);
month = 0;
mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
%
% Initial conditions:
%
if exist(strcat(p.pathN0,'.mat'),'file')
    load(p.pathN0, 'N');
    u(:, ixN) = gridToMatrix(N, [], p.pathBoxes, p.pathGrid);
else
    u(:, ixN) = p.u0(p.idxN)*ones(nb,1);
end
u(:, ixDOC) = zeros(nb,1) + p.u0(ixDOC);
if bSilicate
    if exist(strcat(p.pathSi0,'.mat'),'file')
        load(p.pathSi0, 'Si');
        u(:, ixSi) = gridToMatrix(Si, [], p.pathBoxes, p.pathGrid);
    else
        u(:, ixSi) = zeros(nb,1) + p.u0(ixSi);
    end
end
u(:, ixB) = ones(nb,1)*p.u0(ixB);

if ~isempty(sim)
    if options.bVerbose
        disp('Starting from previous simulation.');
    end
    u(:,ixN) = gridToMatrix(squeeze(double(sim.N(end,:,:,:))),[],p.pathBoxes, p.pathGrid);
    u(:, ixDOC) = gridToMatrix(squeeze(double(sim.DOC(end,:,:,:))),[],p.pathBoxes, p.pathGrid);
    if bSilicate
        u(:, ixSi) = gridToMatrix(squeeze(double(sim.Si(end,:,:,:))),[],p.pathBoxes, p.pathGrid);
    end
    if ~options.bNinit
        for i = 1:p.n -p.idxB+1
            u(:, ixB(i)) = gridToMatrix(squeeze(double(squeeze(sim.B(end,:,:,:,i)))),[],p.pathBoxes, p.pathGrid);
        end
    end
end

sim = load(p.pathGrid,'x','y','z','dznom','bathy'); % Get grid
%
% Load temperature:
%
load(p.pathTemp, 'Tbc');
Tmat = zeros(nb,12);
for i = 1:12
    Tmat(:,i) = gridToMatrix(Tbc(:,:,:,i), [], p.pathBoxes, p.pathGrid);
end
%
% Load Light:
%
if p.bUse_parday_light
    if exist(p.pathPARday,'file')
        load(p.pathPARday);
        for i = 1:length(sim.z)
            parday(:,:,i,:) = parday(:,:,1,:);
        end
        parday = gridToMatrix(parday, [], p.pathBoxes, p.pathGrid);
    else
        error('PARday file does not exist. Set p.bUse_parday_light = false');
    end
end
L0 = zeros(nb,365/p.dtTransport );
for i = 1:365/p.dtTransport
    if p.bUse_parday_light
        % parday has its own time resolution (half days for MITgcm), so pick
        % the column that covers day i*dtTransport:
        iPar = ceil(i*p.dtTransport*size(parday,2)/365);
        L0(:,i) = 1e6*parday(:,iPar)/(24*60*60).*exp(-p.kw*Zbox);
    else
        % Calculate light. i*p.dtTransport is the day of the year; it was written
        % as i/2, which is only the same thing when dtTransport is half a day:
        L0(:,i) = p.EinConv*p.PARfrac*daily_insolation(0,Ybox,i*p.dtTransport,1).*exp(-p.kw*Zbox);
    end
end
%
% Calculate sinking matrices. Uses an explicit first-order upwind scheme:
%
% Get sinking velocities from libNUMmodel:
p.velocity = 0*p.m;
p.velocity = calllib(loadNUMmodelLibrary(), 'f_getsinking', p.velocity);
% Find indices of groups with sinking
idxSinking = find(p.velocity ~= 0);

if ~isempty(idxSinking)
    % Find the indices into the grid
    xx = matrixToGrid((1:nb)', [], p.pathBoxes, p.pathGrid);
    %
    % Every box loses a fraction flx of its mass to the box below and gains the
    % same fraction of the mass of the box above. flx depends only on how deep
    % the box is in its own water column, so the only thing needed per box is
    % that depth and the index of the box above it.
    %
    % The matrices are assembled in one call to sparse() at the end. Writing the
    % elements one at a time into a sparse matrix rebuilds the whole thing on
    % every write, which for the ECCO grid takes minutes rather than a moment.
    %
    [nx, ny, nz] = size(xx);
    X = reshape(xx, nx*ny, nz);
    bValid = ~isnan(X); % The wet boxes
    kDepth = cumsum(bValid, 2); % Depth of each box counted down its own column
    %
    % The box above, which is the previous wet one in the same column:
    %
    iAbove = NaN(size(X));
    iLast = NaN(size(X,1),1);
    for k = 1:nz
        iAbove(:,k) = iLast;
        iLast(bValid(:,k)) = X(bValid(:,k), k);
    end
    bBottom = bValid & (kDepth == kDepth(:,end)); % Deepest box of each column

    ixBox = X(bValid);
    kBox = kDepth(bValid);
    ixAbove = iAbove(bValid);
    bBottomBox = bBottom(bValid);
    bHasAbove = kBox > 1;

    Asink = cell(1,length(idxSinking));
    for l = 1:length(idxSinking)
        flx = min(1, p.velocity(idxSinking(l))*p.dtTransport./sim.dznom(kBox));
        % Loss of mass ...
        vDiagonal = 1 - flx;
        if p.BC_POMclosed
            vDiagonal(bBottomBox) = 1; % Closed BC; no loss of mass at the bottom
        end
        % ... and gain from above
        Asink{l} = sparse( [ixBox; ixBox(bHasAbove)], ...
                           [ixBox; ixAbove(bHasAbove)], ...
                           [vDiagonal; flx(bHasAbove)], nb, nb );
    end
end
%%
% Find indices of bottom cells for bottom boundary condition:
%
xx = matrixToGrid((1:nb)', [], p.pathBoxes, p.pathGrid); % Find the indices into the grid
ixBottom = []; % Indices to all the bottom cells
dzBottom = []; % The width of all the bottom cells.
for i = 1:size(xx,1)
    for j = 1:size(xx,2)
        ixZ = isnan( xx(i,j,:) );
        if ~ixZ(1)
            idx = find(ixZ==0,1,'last'); % Find the last cell
            ixBottom = [ixBottom xx(i,j,idx)];
            dzBottom = [dzBottom sim.dznom(idx)];
        end
    end
end
% Set BCvalue:
BCvalue = p.BCvalue;
if size(BCvalue,1)==1
    BCvalue = ones(length(ixBottom),1)*BCvalue;
end
% If BCvalue == -1 then use the bottom value from the initial conditions:
for i = 1:size(BCvalue,2)
    if BCvalue(1,i)==-1
        BCvalue(:,i) = u(ixBottom,i)';
    end
end
%%
% Water columns for removing negative concentrations after transport
% (function_clip_negative). All boxes are listed column by column, each
% column from the surface down:
%
bClip = isfield(p,'bClipNegative') && p.bClipNegative;
ixClip = []; % The state variables that are clipped
if isfield(p,'bClipNutrients') && p.bClipNutrients
    ixClip = 1:p.nNutrients;
end
if ~isfield(p,'bClipBiomass') || p.bClipBiomass
    ixClip = [ixClip, ixB];
end
%
% The unicellular groups whose deficits may also be taken from neighbouring
% size classes (first and last column within u(:,ixClip)):
%
sizeGroup = zeros(length(ixClip), 2);
if isfield(p,'bClipSizeClasses') && p.bClipSizeClasses
    for iGroup = find(ismember(p.typeGroups, [1 3 4 5])) % Generalists and diatoms
        ix = find(ismember(ixClip, p.ixStart(iGroup):p.ixEnd(iGroup)));
        sizeGroup(ix,:) = repmat([ix(1) ix(end)], length(ix), 1);
    end
end
if bClip
    [nx, ny, nz] = size(xx);
    Xt = reshape(xx, nx*ny, nz)'; % One column per water column
    bWet = ~isnan(Xt);
    clipCols.order = Xt(bWet);
    [~, c] = find(bWet); % Water column of each box
    iCol = cumsum([1; diff(c) ~= 0]);
    firstPos = find([true; diff(c) ~= 0]);
    lastPos = [firstPos(2:end) - 1; numel(c)];
    clipCols.pos = zeros(nb,1);
    clipCols.pos(clipCols.order) = 1:nb;
    clipCols.first = zeros(nb,1);
    clipCols.first(clipCols.order) = firstPos(iCol);
    clipCols.last = zeros(nb,1);
    clipCols.last(clipCols.order) = lastPos(iCol);
end
%%
% Matrices for saving the solution:
%
iSave = 0;
nSave = floor(p.tEnd/p.tSave) + sign(mod(p.tEnd,p.tSave));
sim.N = single(zeros(nSave,length(sim.x), length(sim.y), length(sim.z)));
if bSilicate
    sim.Si = sim.N;
end
sim.DOC = sim.N;
sim.B = single(zeros(nSave, length(sim.x), length(sim.y), length(sim.z), p.n-p.idxB+1));
sim.L = sim.N;
sim.T = sim.N;
tSave = [];
%
% N budget (gN over the whole ocean, double precision). The losses and the
% bottom-BC input are measured as the change in total N over the respective
% sub-steps, accumulated between saves. Whatever is not captured by these is
% the non-conservation of the transport step.
%
S = inputRead;
rhoCN = S.general.rhoCN;
gridDV = load(p.pathGrid,'dv');
dvBox = reshape(gridToMatrix(gridDV.dv, [], p.pathBoxes, p.pathGrid), [], 1); % m3
fNtot = @(uu) ( reshape(uu(:,ixN),1,[]) + reshape(sum(uu(:,ixB),2),1,[])/rhoCN ) * dvBox/1000; % gN
sim.Ntot = zeros(nSave,1);
sim.Nloss = zeros(nSave,1);     % N leaving the ocean by sinking and (without POM) HTL/POM export
sim.NlossHTL = zeros(nSave,1);  % The HTL/POM-export part of Nloss (zero when POM is present)
sim.Nprod = zeros(nSave,1);     % N entering through the bottom BC
sim.Nclip = zeros(nSave,1);     % N added when removing negative concentrations after transport
sim.nutrientsNegative = zeros(nSave,p.nNutrients); % Mass of negative nutrients (gN, gC, gSi), when they are not clipped
sim.Bnegative = zeros(nSave,p.nGroups); % Mass of negative biomass in each group (gC), when it is not clipped
sim.fracBnegative = zeros(nSave,1);     % Fraction of the biomass entries that are negative
NlossBio = 0; NlossSink = 0; Nbc = 0; Nclip = 0;
clipStats = zeros(4, p.n); % Totals from function_clip_negative over the run
%
% Matrices for annual averages:
%
if options.bCalcAnnualAverages
    ProdNet = zeros( simtime, nb,1);
    ProdHTL = zeros( simtime, nb,1);
    mHTL = zeros( simtime, nb,1);
    BHTL = zeros( nb, p.n );
    ProdGrossAnnual = zeros( nb,1 );
    ProdNetAnnual = zeros( nb,1 );
    ProdHTLAnnual = zeros( nb,1 );
    BpicoAnnualMean = zeros( nb,1 );
    BnanoAnnualMean = zeros( nb,1 );
    BmicroAnnualMean = zeros( nb,1 );

    [mortHTL, bQuadratic] = getMortHTL(p);
end

% ---------------------------------------
% Run transport matrix simulation
% ---------------------------------------

if options.bVerbose
    disp('Starting simulation')
end
sLibname = loadNUMmodelLibrary();
%
% The biology is run on the workers if a parallel pool is open. A setup made
% without bParallel=true is only made in this session, and the workers then
% keep simulating whatever setup they had before, so check that they agree:
%
if ~isempty(gcp('nocreate'))
    checkWorkerSetup(p);
end
dtTransport = p.dtTransport;
n = p.n;

tic
for i=1:simtime
    %
    % Test for time to change monthly transport matrix
    %
    if ismember(mod(i,365/p.dtTransport), 1+cumsum(mon)/p.dtTransport)
        % Load TM
        load(strcat(p.pathMatrix, sprintf('%02i.mat',month+1)));
        %disp(strcat(p.pathMatrix, sprintf('%02i.mat',month+1)));

        % With bClip the raw matrices are used and the negative values
        % they create are removed after each transport step:
        if ~bClip
            Aexp = function_convert_TM_positive(Aexp, dvBox);
            Aimp = function_convert_TM_positive(Aimp, dvBox);
        end

        % Preparing for timestepping. 43200s.
        load(p.pathGrid,'deltaT')
        Aexp = Ix + (p.dtTransport*24*60*60)*Aexp;
        Aimp = Aimp^(p.dtTransport*24*60*60/deltaT);

        % Set monthly mean temperature
        T = Tmat(:,month+1);

        month = mod(month + 1, 12);
    end
    %
    % Run Euler time step for dtTransport days:
    %
    L = L0(:,mod(i,365/p.dtTransport)+1);
    dt = p.dt;

    Nbefore = fNtot(u);

    if ~isempty(gcp('nocreate'))
        if options.bCalcAnnualAverages && i > simtime - 365/p.dtTransport
            %
            % Integrate and calculate functions (only last year):
            %
            nNutrients = p.nNutrients;
            parfor k = 1:nb
                ProdGross1 = 0;
                ProdNet1 = 0;
                ProdHTL1 = 0;
                ProdBact1 = 0;
                eHTL1 = 0;
                Bpico1 = 0;
                Bnano1 = 0;
                Bmicro1 = 0;
                mHTL1 = 0;

                [u(k,:), ProdGross1, ProdNet1,ProdHTL1,ProdBact1, eHTL1,Bpico1,Bnano1,Bmicro1,mHTL1] = ...
                    calllib(sLibname, 'f_simulateeulerfunctions', ...
                    u(k,:), L(k), T(k), dtTransport, dt, ...
                    ProdGross1, ProdNet1,ProdHTL1,ProdBact1, eHTL1,Bpico1,Bnano1,Bmicro1,mHTL1);
               
                ProdNet(i,k) = ProdNet1;
                ProdHTL(i,k) = ProdHTL1;
                mHTL(i,k) = mHTL1;
                BHTL(k,:) = BHTL(k,:) + u(k,:).^(1+bQuadratic) .* [zeros(1,nNutrients), mortHTL];

                ProdGrossAnnual(k) = ProdGrossAnnual(k) + ProdGross1;
                ProdNetAnnual(k) = ProdNetAnnual(k) + ProdNet1;
                ProdHTLAnnual(k) = ProdHTLAnnual(k) + ProdHTL1;
                BpicoAnnualMean(k) = BpicoAnnualMean(k) + Bpico1;
                BnanoAnnualMean(k) = BnanoAnnualMean(k) + Bnano1;
                BmicroAnnualMean(k) = BmicroAnnualMean(k) + Bmicro1;
            end
        else
            %
            % Just integrate ODEs:
            %
            parfor k = 1:nb
                u(k,:) = calllib(sLibname, 'f_simulateeuler', ...
                    u(k,:), L(k), T(k), dtTransport, dt);
            end
        end
    else
        for k = 1:nb
            u(k,:) = calllib(sLibname, 'f_simulateeuler', ...
                u(k,:),L(k), T(k), dtTransport, dt);
        end
    end

    Nafter = fNtot(u);
    NlossBio = NlossBio + (Nbefore - Nafter);
    Nbefore = Nafter;
    %
    % Sinking:
    %
    if ~isempty(idxSinking)
        for l = 1:length(idxSinking)
            u(:,idxSinking(l)) = Asink{l}*u(:,idxSinking(l));
        end
    end
    Nafter = fNtot(u);
    NlossSink = NlossSink + (Nbefore - Nafter);
    Nbefore = Nafter;
    %
    % Bottom BC for nutrient fields. The boundary allows diffusion from the
    % bottom into the cell. The diffusivity is controlled by p.BCdiffusion
    % and the bottom value by p.BCvalue. These parameters are set in
    % parametersGlobal:
    %
    for k = 1:p.nNutrients
        u(ixBottom, k) = u(ixBottom, k) +  p.dtTransport* ...
            p.BCmixing(k)./dzBottom'.*(BCvalue(:,k)-u(ixBottom,k));
    end
    Nbc = Nbc + (fNtot(u) - Nbefore);
    %
    % Transport
    %
    if p.bTransport
        u =  Aimp*(Aexp*u);
        if bClip && ~isempty(ixClip)
            st = zeros(4, p.n);
            [u(:,ixClip), st(:,ixClip)] = function_clip_negative(u(:,ixClip), dvBox, clipCols, sizeGroup);
            clipStats = clipStats + st;
            Nclip = Nclip + (st(2,ixN) + sum(st(2,ixB))/rhoCN)/1000; % gN
        end
    end
    %
    % Save timeseries in grid format
    %
    if ((floor(i*(p.dtTransport/p.tSave)) > floor((i-1)*(p.dtTransport/p.tSave))) || (i==simtime))
        if options.bVerbose
            fprintf('t = %u days',floor(i*p.dtTransport))
        end

        if any(isnan(u))
            warning('NaNs after running current time step');
            keyboard
        end

        iSave = iSave + 1;
        sim.N(iSave,:,:,:) = single(matrixToGrid(u(:,ixN), [], p.pathBoxes, p.pathGrid));
        sim.DOC(iSave,:,:,:) = single(matrixToGrid(u(:,ixDOC), [], p.pathBoxes, p.pathGrid));
        if bSilicate
            sim.Si(iSave,:,:,:) = single(matrixToGrid(u(:,ixSi), [], p.pathBoxes, p.pathGrid));
        end
        for j = 1:p.n-p.idxB+1
            sim.B(iSave,:,:,:,j) = single(matrixToGrid(u(:,ixB(j)), [], p.pathBoxes, p.pathGrid));
        end
        sim.L(iSave,:,:,:) = single(matrixToGrid(L, [], p.pathBoxes, p.pathGrid));
        sim.T(iSave,:,:,:) = single(matrixToGrid(T, [], p.pathBoxes, p.pathGrid));
        tSave = [tSave, i*p.dtTransport];
        sim.Ntot(iSave) = fNtot(u);
        sim.Nloss(iSave) = NlossBio + NlossSink;
        sim.NlossHTL(iSave) = NlossBio;
        sim.Nprod(iSave) = Nbc;
        sim.Nclip(iSave) = Nclip;
        sim.nutrientsNegative(iSave,:) = -dvBox'*min(0, u(:,1:p.nNutrients))/1000; % gN, gC, gSi
        for iGroup = 1:p.nGroups
            sim.Bnegative(iSave,iGroup) = -sum(dvBox'*min(0, u(:,p.ixStart(iGroup):p.ixEnd(iGroup))))/1000; % gC
        end
        sim.fracBnegative(iSave) = nnz(u(:,ixB) < 0)/numel(u(:,ixB));
        NlossBio = 0; NlossSink = 0; Nbc = 0; Nclip = 0;
        if options.bVerbose
            fprintf('.\n');
        end
    end
end
time = toc;
if options.bVerbose
    fprintf('Solving time: %2u:%02u:%02u\n', ...
        [floor(time/3600), mod(floor(time/60),60), floor(mod(time,60))]);
end
sim.clipStats = clipStats;
if bClip && options.bVerbose
    %
    % Mass gained by each pool because negative concentrations could not
    % be covered from within the water column. Masses are in mg (u*dv is
    % ugX/l * m3); 1 Tg = 1e15 mg and 1 Pg = 1e18 mg. The percentage is
    % relative to the pool at the end of the run.
    %
    massFinal = dvBox'*u;
    years = p.tEnd/365;
    fprintf('Mass gained from negative concentrations that could not be filled within the water column\n(or, for unicellular groups, from neighbouring size classes):\n');
    fprintf('  %-28s %14s %14s\n', 'Pool', 'Gain per year', '% of pool/yr');
    fPrint = @(name, added, pool, scale, unit) fprintf('  %-28s %9.2e %-4s %14.2e\n', ...
        name, added/years*scale, unit, 100*added/years/pool);
    if ismember(ixN, ixClip)
        fPrint('N', clipStats(2,ixN), massFinal(ixN), 1e-15, 'TgN');
        fPrint('DOC', clipStats(2,ixDOC), massFinal(ixDOC), 1e-18, 'PgC');
        if bSilicate
            fPrint('Si', clipStats(2,ixSi), massFinal(ixSi), 1e-15, 'TgSi');
        end
    end
    bClipB = ismember(ixB(1), ixClip);
    if bClipB
        for iGroup = 1:p.nGroups
            ix = p.ixStart(iGroup):p.ixEnd(iGroup);
            fPrint(p.nameGroup{iGroup}, sum(clipStats(2,ix)), sum(massFinal(ix)), 1e-18, 'PgC');
        end
    end
    NtotFinal = fNtot(u); % gN
    fPrint('Total N (nutrient + biomass)', sum(sim.Nclip), NtotFinal, 1e-12, 'TgN');
    if ~ismember(ixN, ixClip)
        fprintf('Negative nutrients at the end of the run (not clipped):\n');
        fprintf('  %-28s %9.2e TgN  %6.2e %% of pool\n', 'N', ...
            sim.nutrientsNegative(end,ixN)*1e-12, 100*sim.nutrientsNegative(end,ixN)/(massFinal(ixN)/1000));
        fprintf('  %-28s %9.2e PgC  %6.2e %% of pool\n', 'DOC', ...
            sim.nutrientsNegative(end,ixDOC)*1e-15, 100*sim.nutrientsNegative(end,ixDOC)/(massFinal(ixDOC)/1000));
        if bSilicate
            fprintf('  %-28s %9.2e TgSi %6.2e %% of pool\n', 'Si', ...
                sim.nutrientsNegative(end,ixSi)*1e-12, 100*sim.nutrientsNegative(end,ixSi)/(massFinal(ixSi)/1000));
        end
    end
    if ~bClipB
        fprintf('Negative biomass at the end of the run (not clipped; %.2g%% of the entries are negative):\n', ...
            100*sim.fracBnegative(end));
        for iGroup = 1:p.nGroups
            ix = p.ixStart(iGroup):p.ixEnd(iGroup);
            fprintf('  %-28s %9.2e PgC  %6.2e %% of pool\n', p.nameGroup{iGroup}, ...
                sim.Bnegative(end,iGroup)*1e-15, 100*sim.Bnegative(end,iGroup)/(sum(massFinal(ix))/1000));
        end
    end
end

% ---------------------------------------
% Put results into sim structure:
% ---------------------------------------
sim.t = tSave; % days where solution was saved
sim.p = p;
sim.B(sim.B<0) = 0.;
sim.N(sim.N<0) = 0.;
sim.DOC(sim.DOC<0) = 0.;
if bSilicate
    sim.Si(sim.Si<0) = 0.;
end

if options.bCalcAnnualAverages
    sim.ProdNet = zeros(nSave,length(sim.x), length(sim.y));
    sim.ProdHTL = zeros(nSave,length(sim.x), length(sim.y));
    for i = 1:nSave
        ix = floor((simtime/nSave)*(i-1)+1) : (simtime/nSave*i);
        sim.ProdNet(i,:,:) = integrate_over_depth(squeeze(mean(ProdNet(ix,:),1))');
        sim.ProdHTL(i,:,:) = integrate_over_depth(squeeze(mean(ProdHTL(ix,:),1))');
        tmp = single(matrixToGrid( squeeze(mHTL(i,:))', [], p.pathBoxes, p.pathGrid));
        sim.mHTL(i,:,:) = tmp(:,:,1);
    end
    % Average mHTL
    mHTLdepthsum = zeros(length(sim.x), length(sim.y));
    BHTLdepthsum = zeros(length(sim.x), length(sim.y),size(BHTL,2));
    for i = p.idxB:size(BHTL,2)
        BHTLdepthsum(:,:,i) = integrate_over_depth( squeeze(BHTL(:,i)));
        mHTLdepthsum = mHTLdepthsum + log(p.m(i))*BHTLdepthsum(:,:,i);
    end
    sim.mHTLAnnualMean = exp( mHTLdepthsum ./ sum(BHTLdepthsum,3) );
    sim.BHTL = BHTLdepthsum;
    
    % Average productions:
    sim.ProdGrossAnnual(1,:,:) = integrate_over_depth(ProdGrossAnnual);
    sim.ProdNetAnnual(1,:,:) = integrate_over_depth(ProdNetAnnual);
    sim.ProdHTLAnnual(1,:,:) = integrate_over_depth(ProdHTLAnnual);
    % Average biomasses:
    sim.BpicoAnnualMean(1,:,:) = integrate_over_depth(BpicoAnnualMean);
    sim.BnanoAnnualMean(1,:,:) = integrate_over_depth(BnanoAnnualMean);
    sim.BmicroAnnualMean(1,:,:) = integrate_over_depth(BmicroAnnualMean);
end

%%%%%%%%%%%%%%%%%%%%%%
    function area = integrate_over_depth(conc)
        tmp = single(matrixToGrid(conc, [], p.pathBoxes, p.pathGrid)) / (365/p.dtTransport);
        tmp( isnan(tmp) ) = 0;
        area = zeros(length(sim.x), length(sim.y));
        for iz = 1:length(sim.dznom)
            area = area + squeeze(tmp(:,:,iz))*sim.dznom(iz);
        end
    end

end
%
% Checks that all workers of the parallel pool have the same model setup as
% this session, by comparing the derivatives of a test state.
%
function checkWorkerSetup(p)
uTest = reshape(double(p.u0), 1, []);
[~, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', uTest, 100, 15, 0.0, 0*uTest);

f = parfevalOnAll(@workerDerivatives, 1, uTest);
dudtWorkers = fetchOutputs(f);

if ~all(all(abs(dudtWorkers - dudt) <= 1e-12*max(abs(dudt))))
    error('NUMmodel:workerSetup', ['The workers of the parallel pool have a different model setup than this session.\n' ...
        'Make the setup with bParallel=true (e.g. setupNUMmodel(bParallel=true)), ' ...
        'or close the pool with delete(gcp).']);
end
end
%
% Derivatives of the test state on a worker (NaN if its state has another size):
%
function dudt = workerDerivatives(uTest)
sLibname = loadNUMmodelLibrary();
[n, ~] = calllib(sLibname, 'f_getnumgrid', int32(0), int32(0));
if n ~= length(uTest)
    dudt = NaN(size(uTest));
    return
end
[~, dudt] = calllib(sLibname, 'f_calcderivatives', uTest, 100, 15, 0.0, 0*uTest);
end
%
% Removes negative concentrations after a transport step with the raw
% transport matrices. The mass missing in a negative box is taken from the
% nearest boxes above and below in the same water column: first from the
% boxes one step away, then two steps away, etc. The boxes at a given
% distance give in proportion to their mass. If the whole water column
% cannot cover the deficit, the remainder is added as new mass (the box is
% left at zero). No mass is moved between water columns, and each state
% variable is treated separately.
%
% Optionally, a deficit that the water column cannot cover is then taken
% from the neighbouring size classes of the same group in the same box:
% first one size class away, then two, etc., in proportion to their mass.
% That conserves both carbon and nitrogen, and is only meaningful for
% unicellular groups.
%
% Within each water column the negative boxes are corrected from the
% surface down. Positive boxes are only ever scaled down, so no box ends
% up negative.
%
% In:
%  u    - state (nb x n concentrations)
%  dv   - volume of each box (nb x 1, m3)
%  cols - the water columns (made in the main function), with fields:
%         order: all boxes, column by column, each from the surface down
%         pos:   position of each box in order
%         first, last: position in order of the top and bottom box of
%                each box's water column
%  sizeGroup - (optional) n x 2: the first and last column in u of the
%         group that each state variable belongs to, if the deficit may be
%         taken from its neighbouring size classes; [0 0] if not.
%
% Out:
%  u  - state without negative values
%  st - 4 x n statistics for each state variable (concentration * m3):
%       st(1,:): negative mass removed
%       st(2,:): mass added because the water column could not cover it
%       st(3,:): mass taken from boxes more than one step away
%       st(4,:): mass taken from neighbouring size classes
%
function [u, st] = function_clip_negative(u, dv, cols, sizeGroup)

if nargin < 4
    sizeGroup = zeros(size(u,2), 2);
end
st = zeros(4, size(u,2));
[ib, jt] = find(u < 0);
if isempty(ib)
    return
end
%
% Sort the negative entries by state variable and from the surface down:
%
[~, s] = sortrows([jt, cols.pos(ib)]);
ib = ib(s);
jt = jt(s);

for m = 1:numel(ib)
    b = ib(m);
    j = jt(m);
    d = -u(b,j)*dv(b); % The deficit
    st(1,j) = st(1,j) + d;
    u(b,j) = 0;

    q = cols.pos(b);
    q0 = cols.first(b);
    q1 = cols.last(b);
    r = 1;
    while d > 0 && (q-r >= q0 || q+r <= q1)
        % The positive boxes r steps above and below:
        bAbove = 0;
        bBelow = 0;
        avail = 0;
        if q-r >= q0
            bAbove = cols.order(q-r);
            if u(bAbove,j) > 0
                avail = avail + u(bAbove,j)*dv(bAbove);
            else
                bAbove = 0;
            end
        end
        if q+r <= q1
            bBelow = cols.order(q+r);
            if u(bBelow,j) > 0
                avail = avail + u(bBelow,j)*dv(bBelow);
            else
                bBelow = 0;
            end
        end
        % Take the deficit in proportion to their mass:
        if avail > 0
            take = min(d, avail);
            f = 1 - take/avail;
            if bAbove > 0
                u(bAbove,j) = u(bAbove,j)*f;
            end
            if bBelow > 0
                u(bBelow,j) = u(bBelow,j)*f;
            end
            if r > 1
                st(3,j) = st(3,j) + take;
            end
            d = d - take;
        end
        r = r + 1;
    end
    %
    % Take what is missing from the neighbouring size classes in the box:
    %
    j0 = sizeGroup(j,1);
    j1 = sizeGroup(j,2);
    r = 1;
    while d > 0 && j0 > 0 && (j-r >= j0 || j+r <= j1)
        avail = 0;
        jj = [j-r, j+r];
        jj = jj(jj >= j0 & jj <= j1);
        jj = jj(u(b,jj) > 0);
        if ~isempty(jj)
            avail = sum(u(b,jj))*dv(b);
        end
        if avail > 0
            take = min(d, avail);
            u(b,jj) = u(b,jj)*(1 - take/avail);
            st(4,j) = st(4,j) + take;
            d = d - take;
        end
        r = r + 1;
    end
    st(2,j) = st(2,j) + d; % What could not be covered is added
end
end
