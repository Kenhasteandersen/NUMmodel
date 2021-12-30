%
% Global run using transport matrices
% 
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% (choose MITgcm_2.8deg), and be put into the location 'NUMmodel/TMs'
%
% Input:
%  p: parameter structure from parametersGlobal
%  sim: (optional) simulation to use for initial conditions
%  bCalcAnnualAverages: increases the simulation time by a factor 2-3
%
% Output:
%  sim: structure with simulation results
%
function sim = simulateGlobal(p, sim, bCalcAnnualAverages)

arguments
    p struct;
    sim struct = [];
    bCalcAnnualAverages = false; % Whether to calculate annual averages
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
% Check that files exist:
%
path = fileparts(mfilename('fullpath'));
addpath(strcat(path,'/Transport matrix'));

if ~exist(p.pathBoxes,'file')
    error( sprintf('Error: Cannot find transport matrix file: %s',...
        p.pathBoxes));
end
% ---------------------------------
% Load library:
% ---------------------------------
% if p.bParallel
%     if isempty(gcp('nocreate'))
%         parpool('AttachedFiles',...
%             {'../Fortran/NUMmodel_matlab.so',...
%              '../Fortran/NUMmodel_wrap_colmajor4matlab.h'});
%     end
%     %
%     % Set parameters:
%     %
%     h = gcp('nocreate');
%     poolsize = h.NumWorkers;
%     parfor i=1:poolsize
%         loadNUMmodelLibrary();
%         %calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(p.mAdult)), p.mAdult);
%         calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly',int32(10));
%     end
% else
%     loadNUMmodelLibrary();
%     %calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(p.mAdult)), p.mAdult);
%     calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly',int32(10));
% end
% ---------------------------------------
% Initialize run:
% ---------------------------------------
simtime = p.tEnd*2; %simulation time in half days
% Load grid data:
%load(p.pathGrid);
%load(p.pathConfigData);
load(p.pathBoxes, 'nb', 'Ybox', 'Zbox');

% Preparing timestepping
Ix = speye(nb,nb);
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
  load 'Transport Matrix/parday';
end
L0 = zeros(nb,730);
for i = 1:730
    if p.bUse_parday_light
        L0(:,i) = 1e6*parday(:,i)/(24*60*60).*exp(-p.kw*Zbox);
    else
        % Calculate light:
        L0(:,i) = p.EinConv*p.PARfrac*daily_insolation(0,Ybox,i/2,1).*exp(-p.kw*Zbox);
    end
end
%
% Matrices for saving the solution:
%
iSave = 0;
nSave = floor(p.tEnd/p.tSave) + sign(mod(p.tEnd,p.tSave));
sim = load(p.pathGrid,'x','y','z','dznom','bathy');
sim.N = single(zeros(length(sim.x), length(sim.y), length(sim.z),nSave));
if bSilicate
    sim.Si = sim.N;
end
sim.DOC = sim.N;
sim.B = single(zeros(length(sim.x), length(sim.y), length(sim.z), p.n-p.idxB+1, nSave));
sim.L = sim.N;
sim.T = sim.N;
tSave = [];
%
% Matrices for annual averages:
%
if bCalcAnnualAverages
   sim.ProdGrossAnnual = zeros( nb,1 );
   sim.ProdNetAnnual = zeros( nb,1 );
   sim.ProdHTLAnnual = zeros( nb,1 );
   sim.BpicoAnnualMean = zeros( nb,1 );
   sim.BnanoAnnualMean = zeros( nb,1 );
   sim.BmicroAnnualMean = zeros( nb,1 );
end

% ---------------------------------------
% Run transport matrix simulation
% ---------------------------------------
disp('Starting simulation')
sLibname = loadNUMmodelLibrary();

tic
for i=1:simtime
    %
    % Test for time to change monthly transport matrix
    %
    if ismember(mod(i,730), 1+2*cumsum(mon))
        % Load TM
        load(strcat(p.pathMatrix, sprintf('%02i.mat',month+1)));
        %disp(strcat(p.pathMatrix, sprintf('%02i.mat',month+1)));
        
        Aexp=function_convert_TM_positive(Aexp);
        Aimp=function_convert_TM_positive(Aimp);
        
        % Preparing for timestepping. 43200s.
        load(p.pathGrid,'deltaT')
        Aexp = Ix + (12*60*60)*Aexp;
        Aimp = Aimp^(12*60*60/deltaT);
        
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
        parfor k = 1:nb
            u(k,:) = calllib(sLibname, 'f_simulateeuler', ...
                u(k,:), L(k), T(k), 0.5, dt);
        end
    else
        for k = 1:nb
            u(k,:) = calllib(sLibname, 'f_simulateeuler', ...
                u(k,:),L(k), T(k), 0.5, dt);
            %u(k,1) = u(k,1) + 0.5*(p.u0(1)-u(k,1))*0.5;
            % If we use ode23:
            %[t, utmp] = ode23(@fDerivLibrary, [0 0.5], u(k,:), [], L(k));
            %u(k,:) = utmp(end,:);
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
        %for k = 1:p.n
        %    u(:,k) =  Aimp * (Aexp * u(:,k));
        %end
        u =  Aimp*(Aexp*u);
    end
    %
    % Enforce minimum B concentration
    %
    u(u<p.umin) = p.umin; 
    
    %
    % Save timeseries in grid format
    %
    if ((mod(i/2,p.tSave) < mod((i-1)/2,p.tSave)) || (i==simtime))
        fprintf('t = %u days',floor(i/2))
        iSave = iSave + 1;
        sim.N(:,:,:,iSave) = single(matrixToGrid(u(:,ixN), [], p.pathBoxes, p.pathGrid));
        sim.DOC(:,:,:,iSave) = single(matrixToGrid(u(:,ixDOC), [], p.pathBoxes, p.pathGrid));
        if bSilicate
            sim.Si(:,:,:,iSave) = single(matrixToGrid(u(:,ixSi), [], p.pathBoxes, p.pathGrid));
        end
        for j = 1:p.n-p.idxB+1
            sim.B(:,:,:,j,iSave) = single(matrixToGrid(u(:,ixB(j)), [], p.pathBoxes, p.pathGrid));
        end
        sim.L(:,:,:,iSave) = single(matrixToGrid(L, [], p.pathBoxes, p.pathGrid));
        sim.T(:,:,:,iSave) = single(matrixToGrid(T, [], p.pathBoxes, p.pathGrid));
        tSave = [tSave, i*0.5];
        fprintf('.\n');      
    end
    %
    % Update annual averages:
    %
    if bCalcAnnualAverages
        for k = 1:nb
            [ProdGross1, ProdNet1,ProdHTL1,eHTL,Bpico1,Bnano1,Bmicro1] = ...
                            getFunctions(u(k,:), L(k), T(k));
            sim.ProdGrossAnnual(k) = sim.ProdGrossAnnual(k) + ProdGross1/(p.tEnd*2);            
            sim.ProdNetAnnual(k) = sim.ProdNetAnnual(k) + ProdNet1/(p.tEnd*2);
            sim.ProdHTLAnnual(k) = sim.ProdHTLAnnual(k) + ProdHTL1/(p.tEnd*2);
            sim.BpicoAnnualMean(k) = sim.BpicoAnnualMean(k) + Bpico1/(p.tEnd*2*365);
            sim.BnanoAnnualMean(k) = sim.BnanoAnnualMean(k) + Bnano1/(p.tEnd*2*365);
            sim.BmicroAnnualMean(k) = sim.BmicroAnnualMean(k) + Bmicro1/(p.tEnd*2*365);
        end
    end
    
end
time = toc;
fprintf('Solving time: %2u:%02u:%02u\n', ...
    [floor(time/3600), mod(floor(time/60),60), floor(mod(time,60))]);
% ---------------------------------------
% Put results into sim structure:
% ---------------------------------------
sim.t = tSave; % days where solution was saved
sim.p = p;
sim.Ntot = calcGlobalN(sim);
sim.B(sim.B<0) = 0.;
sim.DOC(sim.DOC<0) = 0.;

if bCalcAnnualAverages
    tmp = single(matrixToGrid(sim.ProdGrossAnnual, [], p.pathBoxes, p.pathGrid));
    sim.ProdGrossAnnual = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.ProdNetAnnual, [], p.pathBoxes, p.pathGrid));
    sim.ProdNetAnnual = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.ProdHTLAnnual, [], p.pathBoxes, p.pathGrid));
    sim.ProdHTLAnnual = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.BpicoAnnualMean, [], p.pathBoxes, p.pathGrid));
    sim.BpicoAnnualMean = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.BnanoAnnualMean, [], p.pathBoxes, p.pathGrid));
    sim.BnanoAnnualMean = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.BmicroAnnualMean, [], p.pathBoxes, p.pathGrid));
    sim.BmicroAnnualMean = squeeze(tmp(:,:,1));
end
end
