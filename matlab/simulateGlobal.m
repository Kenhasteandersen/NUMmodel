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
%  sim: (optional) simulation to use for initial conditions
%  options.bCalcAnnualAverages: increases the simulation time of the last year by a
%                       factor 10
%  options.bVerbose: Write out progress on the terminal
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
% ---------------------------------------
% Initialize run:
% ---------------------------------------
simtime = p.tEnd/p.dtTransport; %simulation time in half days
load(p.pathBoxes, 'nb', 'Ybox', 'Zbox');

% Preparing timestepping
Ix = speye(nb,nb);
month = 0;
mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
%
% Initial conditions:
%
if ~isempty(sim)
    if options.bVerbose
        disp('Starting from previous simulation.');
    end
    u(:,ixN) = gridToMatrix(squeeze(double(sim.N(end,:,:,:))),[],p.pathBoxes, p.pathGrid);
    u(:, ixDOC) = gridToMatrix(squeeze(double(sim.DOC(end,:,:,:))),[],p.pathBoxes, p.pathGrid);
    if bSilicate
        u(:, ixSi) = gridToMatrix(squeeze(double(sim.Si(end,:,:,:))),[],p.pathBoxes, p.pathGrid);
    end
    for i = 1:p.n -p.idxB+1
        u(:, ixB(i)) = gridToMatrix(squeeze(double(squeeze(sim.B(end,:,:,:,i)))),[],p.pathBoxes, p.pathGrid);
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
        if exist(strcat(p.pathSi0,'.mat'),'file')
            load(p.pathSi0, 'Si');
            u(:, ixSi) = gridToMatrix(Si, [], p.pathBoxes, p.pathGrid);
        else
            u(:, ixSi) = zeros(nb,1) + p.u0(ixSi);
        end
    end
    u(:, ixB) = ones(nb,1)*p.u0(ixB);
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
        L0(:,i) = 1e6*parday(:,i)/(24*60*60).*exp(-p.kw*Zbox);
    else
        % Calculate light:
        L0(:,i) = p.EinConv*p.PARfrac*daily_insolation(0,Ybox,i/2,1).*exp(-p.kw*Zbox);
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
    if options.bVerbose
        disp('Allocating sinking matrices')
    end
    % Allocate sinking matrices:
    Asink = {};
    for l = 1:length(idxSinking)
        Asink{l} = sparse(1,1,0,nb,nb,nb*2);
    end
    % Find the indices into the grid
    xx = matrixToGrid((1:nb)', [], p.pathBoxes, p.pathGrid);
    % Run through all latitudes and longitudes:
    for i = 1:size(xx,1)
        for j = 1:size(xx,2)
            % Find the watercolumn indices:
            idxGrid = squeeze(xx(i, j, :));
            idxGrid = idxGrid( ~isnan(idxGrid));
            if ~isempty(idxGrid)
                % Run through all sinking state variables
                for l = 1:length(idxSinking)
                    for k = 1:length(idxGrid)
                        flx = min(1, p.velocity(idxSinking(l))*p.dtTransport./sim.dznom(k));
                        % Loss of mass ...
                        Asink{l}(idxGrid(k),idxGrid(k)) = 1-flx;
                        % Gain from above
                        if (k > 1)
                            Asink{l}(idxGrid(k),idxGrid(k-1)) = flx;
                        end
                    end
                    if p.BC_POMclosed
                        Asink{l}(idxGrid(k),idxGrid(k)) = 1; % Closed BC; no loss of mass at the bottom
                    end
                end
            end
        end
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
% Matrices for annual averages:
%
if options.bCalcAnnualAverages
    ProdGrossAnnual = zeros( nb,1 );
    ProdNetAnnual = zeros( nb,1 );
    ProdHTLAnnual = zeros( nb,1 );
    BpicoAnnualMean = zeros( nb,1 );
    BnanoAnnualMean = zeros( nb,1 );
    BmicroAnnualMean = zeros( nb,1 );
end

% ---------------------------------------
% Run transport matrix simulation
% ---------------------------------------
if options.bVerbose
    disp('Starting simulation')
end
sLibname = loadNUMmodelLibrary();
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

        Aexp = function_convert_TM_positive(Aexp);
        Aimp = function_convert_TM_positive(Aimp);

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

    %N = calc_tot_n(p,u);

    if ~isempty(gcp('nocreate'))
        if options.bCalcAnnualAverages && i > simtime - 365/p.dtTransport
            %
            % Integrate and calculate functions (only last year):
            %
            parfor k = 1:nb
                ProdGross1 = 0;
                ProdNet1 = 0;
                ProdHTL1 = 0;
                ProdBact1 = 0;
                eHTL1 = 0;
                Bpico1 = 0;
                Bnano1 = 0;
                Bmicro1 = 0;

                [u(k,:), ProdGross1, ProdNet1,ProdHTL1,ProdBact1, eHTL1,Bpico1,Bnano1,Bmicro1] = ...
                    calllib(sLibname, 'f_simulateeulerfunctions', ...
                    u(k,:), L(k), T(k), dtTransport, dt, ...
                    ProdGross1, ProdNet1,ProdHTL1,ProdBact1, eHTL1,Bpico1,Bnano1,Bmicro1);

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

    %calc_tot_n(p,u)/N - 1
    %N = calc_tot_n(p,u);
    %
    % Sinking:
    %
    if ~isempty(idxSinking)
        for l = 1:length(idxSinking)
            u(:,idxSinking(l)) = Asink{l}*u(:,idxSinking(l));
        end
    end
    %calc_tot_n(p,u)/N - 1
    %N = calc_tot_n(p,u);
    %
    % Bottom BC for nutrient fields. The boundary allows diffusion from the
    % bottom into the cell. The diffusivity is controlled by p.BCmixing
    % and the bottom value by p.BCvalue. These parameters are set in
    % parametersGlobal:
    %
    for k = 1:p.nNutrients
        u(ixBottom, k) = u(ixBottom, k) +  p.dtTransport* ...
            p.BCmixing(k).*(BCvalue(:,k)-u(ixBottom,k));
    end
    %calc_tot_n(p,u)/N - 1
    %N = calc_tot_n(p,u);
    %
    % Transport
    %
    if p.bTransport
        u =  Aimp*(Aexp*u);
    end
    %calc_tot_n(p,u)/N - 1
    %N = calc_tot_n(p,u);
    %fprintf('---\n')
    %
    % Save timeseries in grid format
    %
    if ((floor(i*(p.dtTransport/p.tSave)) > floor((i-1)*(p.dtTransport/p.tSave))) || (i==simtime))
        if options.bVerbose
            fprintf('t = %u days',floor(i/2))
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

% ---------------------------------------
% Put results into sim structure:
% ---------------------------------------
sim.t = tSave; % days where solution was saved
sim.p = p;
sim.Ntot = calcGlobalN(sim);
sim.B(sim.B<0) = 0.;
sim.DOC(sim.DOC<0) = 0.;

if options.bCalcAnnualAverages
    sim.ProdGrossAnnual(1,:,:) = integrate_over_depth(ProdGrossAnnual);
    sim.ProdNetAnnual(1,:,:) = integrate_over_depth(ProdNetAnnual);
    sim.ProdHTLAnnual(1,:,:) = integrate_over_depth(ProdHTLAnnual);
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

