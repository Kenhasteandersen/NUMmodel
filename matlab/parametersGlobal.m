%
% Sets the parameters for the global model simlations
%
% In:
%  n - number of size groups (default 10)
%  nTMmodel - Which transport matrices to use:
%       1 = MITgcm_2.8 (default)
%       2 = MITgcm_ECCO
%       3 = UVicOSUpicdefault (experimental)
%
% Out:
%  simulation structure
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
    
    p.dt = 0.1; % For Euler time stepping
elseif nTMmodel == 2
    p.TMname = 'MITgcm_ECCO';
    p.pathMatrix = strcat(path,'/../TMs/MITgcm_ECCO/Matrix1/TMs/matrix_nocorrection_');
    p.pathBoxes = strcat(path,'/../TMs/MITgcm_ECCO/Matrix1/Data/boxes.mat');
    p.pathGrid = strcat(path,'/../TMs/MITgcm_ECCO/grid.mat');
    p.pathConfigData = strcat(path,'/../TMs/MITgcm_ECCO/config_data.mat');
    p.pathTemp = strcat(path,'/../TMs/MITgcm_ECCO/BiogeochemData/Theta_bc.mat'); 
    p.pathN0    = strcat(path,'/../TMs/MITgcm_ECCO_N0');
    p.pathInit = strcat(sprintf('Transport matrix/globalInitMITgcm_ECCO_%02i',length(p.u0)));
    
    p.dt = 0.1; % For Euler time stepping
elseif nTMmodel == 3
    % Experimental
     p.TMname = 'UVicOSUpicdefault';
    p.pathMatrix = strcat(path,'/../TMs/UVicOSUpicdefault/Matrix1/TMs/matrix_nocorrection_');
    p.pathBoxes = strcat(path,'/../TMs/UVicOSUpicdefault/Matrix1/Data/boxes.mat');
    p.pathGrid = strcat(path,'/../TMs/UVicOSUpicdefault/grid.mat');
    p.pathConfigData = strcat(path,'/../TMs/UVicOSUpicdefault/config_data.mat');
    p.pathTemp = strcat(path,'/../TMs/UVicOSUpicdefault/BiogeochemData/Theta_bc.mat'); 
    p.pathN0    = strcat(path,'/../TMs/UVicOSUpicdefault_N0');
    p.pathInit = strcat(sprintf('Transport matrix/globalInitUVicOSUpicdefault_%02i',length(p.u0)));
    
    p.dt = 0.1; % For Euler time stepping

end
%
% Numerical parameters:
%
p.tEnd = 365; % In days
p.tSave = 365/12; % How often to save results (monthly)
p.bTransport = true;
p.umin = 1e-5*p.mDelta(3)/p.m(3); % Minimum B concentration
%
% Light environment:
%
p.bUse_parday_light = false; % Using the parday file includes changes in cloud cover
                             % but only works with MITgcm_2.8
% Parameters used to calculate light if not using parday:
p.EinConv = 4.57; % conversion factor from W m^-2 to \mu mol s^-1 m^-2 (Thimijan & Heins 1983)
p.PARfrac = 0.4; % Fraction of light available as PAR. Source unknown
p.kw = 0.05; % m^-1
