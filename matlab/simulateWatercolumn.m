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
%  bCalcAnnualAverages: increases the simulation time by a factor 2-3
%
% Output:
%  sim: structure with simulation results
%
function sim = simulateWatercolumn(p, lat, lon, sim)

arguments
    p struct;
    lat = 0;
    lon = 0;
    sim struct = [];
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

% ---------------------------------------
% Find the indices that corresponds to the selected water column:
% ---------------------------------------
sim = load(p.pathGrid,'x','y','z','dznom','bathy'); % Load grid
load(p.pathBoxes, 'nb', 'Ybox', 'Zbox');
idx = calcGlobalWatercolumn(lat, lon, sim); % Find the indices into matrix
xx = matrixToGrid((1:nb)', [], p.pathBoxes, p.pathGrid); % Find the indices into the grid
idxGrid = squeeze(xx(idx.x, idx.y, idx.z));
nGrid = length(idxGrid);
%
% Check if a file with the water column exists; if not extract from TM and
% save:
%
Ix = speye(nb,nb);
sFile = sprintf('Watercolumn_%s_lat%03i_lon%02i.mat',p.TMname,lat,lon);
if exist(sFile,'file')
    load(sFile);
else
    fprintf('-> Extracting water column from transport matrix')
    %
    % Check that transport matrix files exist:
    %
    path = fileparts(mfilename('fullpath'));
    addpath(strcat(path,'/Transport matrix'));
    
    if ~exist(p.pathBoxes,'file')
        error( sprintf('Error: Cannot find transport matrix file: %s',...
            p.pathBoxes));
    end
    % Load TMs
    for month=1:12
        load(strcat(p.pathMatrix, sprintf('%02i.mat',month)));
        %disp(strcat(p.pathMatrix, sprintf('%02i.mat',month+1)));
        
        AexpM(month,:,:) = full(function_convert_TM_positive(Aexp(idxGrid,idxGrid)));
        AimpM(month,:,:) = full(function_convert_TM_positive(Aimp(idxGrid,idxGrid)));
        
        % Preparing for timestepping. 43200s.
        AexpM(month,:,:) = Ix(idxGrid,idxGrid) + squeeze((12*60*60)*AexpM(month,:,:));
        AimpM(month,:,:) = squeeze(AimpM(month,:,:))^(36);
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
    % Load Light:
    %
    L0 = zeros(nGrid,730);
    for i = 1:730
        L0(:,i) = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid),i/2,1).*exp(-p.kw*Zbox(idxGrid));
    end
    fprintf('\n');
    save(sFile,'AexpM','AimpM','Tmat','L0');
end

% ---------------------------------------
% Initialize run:
% ---------------------------------------
simtime = p.tEnd*2; %simulation time in half days

% Preparing timestepping

month = 0;
mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
%
% Initial conditions:
%
if (nargin==2)
    disp('Starting from previous simulation.');
    u(:,ixN) = gridToMatrix(squeeze(double(sim.N(:,:,:,end))),[],p.pathBoxes, p.pathGrid);
    u(:, ixDOC) = gridToMatrix(squeeze(double(sim.DOC(:,:,:,end))),[],p.pathBoxes, p.pathGrid);
    if bSilicate
        u(:, ixSi) = gridToMatrix(squeeze(double(sim.Si(:,:,:,end))),[],p.pathBoxes, p.pathGrid);
    end
    for i = 1:p.n -p.idxB+1
        u(:, ixB(i)) = gridToMatrix(squeeze(double(squeeze(sim.B(:,:,:,i,end)))),[],p.pathBoxes, p.pathGrid);
    end
else
    if exist(strcat(p.pathN0,'.mat'),'file')
        load(p.pathN0, 'N');
        u(:, ixN) = gridToMatrix(N, [], p.pathBoxes, p.pathGrid);
    else
        u(:, ixN) = 150*ones(nb,1);
    end
    u(:, ixDOC) = zeros(nb,1) + p.u0(ixDOC);
    if bSilicate
        u(:, ixSi) = zeros(nb,1) + p.u0(ixSi);
    end
    u(:, ixB) = ones(nb,1)*p.u0(ixB);
end
u = u(idxGrid,:); % Use only the specific water column
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
for i=1:simtime
    %
    % Test for time to change monthly temperature
    %
    if ismember(mod(i,730), 1+2*cumsum(mon))
        fprintf('t = %u days\n',floor(i/2))
        % Set monthly mean temperature
        T = Tmat(:,month+1);
        month = mod(month + 1, 12);
    end
    %
    % Run Euler time step for half a day:
    %
    L = L0(:,mod(i,365*2)+1);
    dt = p.dt;
    n = p.n;
    if ~isempty(gcp('nocreate'))
        parfor k = 1:nGrid
            u(k,:) = calllib(loadNUMmodelLibrary(), 'f_simulateeuler', ...
                int32(n), u(k,:), L(k), T(k), 0.5, dt);
        end
    else
        for k = 1:nGrid
            u(k,:) = calllib(loadNUMmodelLibrary(), 'f_simulateeuler', ...
                int32(n), u(k,:),L(k), T(k), 0.5, dt);
        end
    end
    
    if any(isnan(u))
        warning('NaNs after running current time step');
        keyboard
    end
    %
    % Transport
    %
    if p.bTransport
        for k = 1:p.n
            u(:,k) =  squeeze(AimpM(month+1,:,:)) * (squeeze(AexpM(month+1,:,:)) * u(:,k));
            % Enforce conservation crudely:
            u(:,k) = u(:,k) ./(squeeze(AimpM(month+1,:,:)) * (squeeze(AexpM(month+1,:,:)) * ones(nGrid,1)));
        end
    end
    %
    % Enforce minimum concentraion
    %
    u(u<p.umin) = p.umin;
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
        tSave = [tSave, i*0.5];
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
