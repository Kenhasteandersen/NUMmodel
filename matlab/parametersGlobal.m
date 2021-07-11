%
% Sets the parameters for the global model simlations
%
% Input:
%  n - number of size groups (default 10)
%  nTMmodel - Which transport matrices to use:
%       1 = MITgcm_2.8 (default)
%       2 = MITgcm_ECCO
%
function p = parametersGlobal(p, nTMmodel)

arguments
    p struct
    nTMmodel {mustBeInteger} = 1;
end

p.nameModel = 'global';

path = fileparts(mfilename('fullpath'));

%
% Set load paths for tranport matrices:
%
if (nargin==1 || nargin==0 || nTMmodel == 1)
    p.TMname = 'MITgcm_2.8';
    p.pathMatrix   = strcat(path,'/../TMs/MITgcm_2.8deg/Matrix5/TMs/matrix_nocorrection_');
    p.pathBoxes     = strcat(path,'/../TMs/MITgcm_2.8deg/Matrix5/Data/boxes.mat');
    p.pathGrid      = strcat(path,'/../TMs/MITgcm_2.8deg/grid.mat');
    p.pathConfigData = strcat(path,'/../TMs/MITgcm_2.8deg/config_data.mat');
    p.pathTemp      = strcat(path,'/../TMs/MITgcm_2.8deg/BiogeochemData/Theta_bc.mat'); 
    p.pathN0        = strcat(path,'/../TMs/MITgcm_N0');
    p.pathInit      = strcat(sprintf('Transport matrix/globalInitMITgcm_%02i',length(p.u0)));
elseif nTMmodel == 2
    p.TMname = 'MITgcm_ECCO';
    p.pathMatrix = strcat(path,'/../TMs/MITgcm_ECCO/Matrix1/TMs/matrix_nocorrection_');
    p.pathBoxes = strcat(path,'/../TMs/MITgcm_ECCO/Matrix1/Data/boxes.mat');
    p.pathGrid = strcat(path,'/../TMs/MITgcm_ECCO/grid.mat');
    p.pathConfigData = strcat(path,'/../TMs/MITgcm_ECCO/config_data.mat');
    p.pathTemp = strcat(path,'/../TMs/MITgcm_ECCO/BiogeochemData/Theta_bc.mat'); 
    p.pathN0    = strcat(path,'/../TMs/MITgcm_ECCO_N0');
    p.pathInit = strcat(sprintf('Transport matrix/globalInitMITgcm_ECCO_%02i',length(p.u0)));
end
%
% Numerical parameters:
%
p.tEnd = 365; % In days
p.tSave = 365/12; % How often to save results (monthly)
p.dt = 0.1; % For Euler time stepping
p.bTransport = true;
p.umin = 1e-5*p.mDelta(3)/p.m(3); % Minimum B concentration
%
% Light environment (??):
%
p.EinConv = 4.57; % conversion factor from W m^-2 to \mu mol s^-1 m^-2 (Thimijan & Heins 1983)
p.PARfrac = 0.4; % ??
p.kw = 0.1; % 0.4 Camila ??
