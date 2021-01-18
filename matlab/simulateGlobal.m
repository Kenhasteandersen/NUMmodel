%
% Global run of size based unicelluar plankton model
% 
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% (choose MITgcm_2.8deg), and be put into the location '../TMs'
%
% ../Cpp/model.cpp must be compiled locally to create model.so
%
% Input:
%  p: parameter structure from parametersGlobal
%  sim: (optional) simulation to use for initial conditions
%
% Output:
%  sim: structure with simulation results
%
function sim = simulateGlobal(p, sim)
ixN = 1;
ixDOC = 2;
ixB = 3:(2+p.nGrid);

%Tbc = [];

disp('Preparing simulation')
% ---------------------------------
% Load library:
% ---------------------------------
if p.bParallel
    if isempty(gcp('nocreate'))
        parpool('AttachedFiles',...
            {'../Fortran/NUMmodel.so',...
            '../Fortran/NUMmodel_wrap_colmajor4matlab.h'});
    end
    %
    % Set parameters:
    %
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        loadNUMmodelLibrary();
        %calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(mAdult)), mAdult);
        calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly');
    end
else
    loadNUMmodelLibrary();
%    calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(p.mAdult)), p.mAdult);
    calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly');
end
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
    u(:,ixN) = gridToMatrix(squeeze(double(sim.N(:,:,:,end))),[],sim.p.pathBoxes, sim.p.pathGrid);
    u(:, ixDOC) = gridToMatrix(squeeze(double(sim.DOC(:,:,:,end))),[],sim.p.pathBoxes, sim.p.pathGrid);
    for i = 1:p.nGrid
        u(:, ixB(i)) = gridToMatrix(squeeze(double(squeeze(sim.B(:,:,:,i,end)))),[],sim.p.pathBoxes, sim.p.pathGrid);
    end
else
    load(p.pathN0, 'N');
    u(:, ixN) = gridToMatrix(N, [], p.pathBoxes, p.pathGrid);
    u(:, ixDOC) = zeros(nb,1) + p.u0(ixDOC);
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
L0 = zeros(nb,730);
for i = 1:730
    L0(:,i) = p.EinConv*p.PARfrac*daily_insolation(0,Ybox,i/2,1).*exp(-p.kw*Zbox);
end
%
% Matrices for saving the solution:
%
iSave = 0;
nSave = floor(p.tEnd/p.tSave) + sign(mod(p.tEnd,p.tSave));
sim = load(p.pathGrid,'x','y','z');
sim.N = single(zeros(length(sim.x), length(sim.y), length(sim.z),nSave));
sim.DOC = sim.N;
sim.B = single(zeros(length(sim.x), length(sim.y), length(sim.z), p.nGrid, nSave));
sim.L = sim.N;
sim.T = sim.N;
tSave = [];
% ---------------------------------------
% Run transport matrix simulation
% ---------------------------------------
disp('Starting simulation')
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
        Aexp = Ix + (12*60*60)*Aexp;
        Aimp = Aimp^(36);
        
        % Set monthly mean temperature
        T = Tmat(:,month+1);
        
        month = mod(month + 1, 12);
    end
    %
    % Run Euler time step for half a day:
    %
    L = L0(:,mod(i,365*2)+1);
    dt = p.dt;
    if p.bParallel
        parfor k = 1:nb
            u(k,:) = calllib(loadNUMmodelLibrary(), 'f_simulateeuler', int32(p.nGrid+2), u(k,:), L(k), 0.5, dt);
        end
    else
        for k = 1:nb
            u(k,:) = calllib(loadNUMmodelLibrary(), 'f_simulateeuler', int32(p.nGrid+2), u(k,:), L(k), 0.5, dt);
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
        % Could perhaps be parfor'ed for large matrices:
        for k = 1:p.nGrid+2
            u(:,k) =  Aimp * (Aexp * u(:,k));
        end
    end
    %
    % Save timeseries in grid format
    %
    if ((mod(i/2,p.tSave) < mod((i-1)/2,p.tSave)) || (i==simtime))
        fprintf('t = %u days',floor(i/2))
        iSave = iSave + 1;
        sim.N(:,:,:,iSave) = single(matrixToGrid(u(:,ixN), [], p.pathBoxes, p.pathGrid));
        sim.DOC(:,:,:,iSave) = single(matrixToGrid(u(:,ixDOC), [], p.pathBoxes, p.pathGrid));
        for j = 1:p.nGrid
            sim.B(:,:,:,j,iSave) = single(matrixToGrid(u(:,ixB(j)), [], p.pathBoxes, p.pathGrid));
        end
        sim.L(:,:,:,iSave) = single(matrixToGrid(L, [], p.pathBoxes, p.pathGrid));
        sim.T(:,:,:,iSave) = single(matrixToGrid(T, [], p.pathBoxes, p.pathGrid));
        tSave = [tSave, i*0.5];
        fprintf('.\n');      
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

%
% Function to assemble derivative for chemostat:
%
%    function dudt = fDerivLibrary(t,u,L)
%        dudt = 0*u';
%        [u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', length(u), u, L, 0.0, dudt);
%        dudt = dudt';
%    end

end
