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

    function check(sFilename)
       if ~exist(sFilename, 'file')
           error('  Did not find the file:\n  p.TMname%s.\n  Check that transport matrices are downloaded and placed in ../TMs/%s.\nDownload TMs from \n  http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/',...
               sFilename,p.TMname);
       end
    end

p.nameModel = 'global';

path = fileparts(mfilename('fullpath'));

%
% Set load paths for tranport matrices:
%
if (nargin==1 || nargin==0 || nTMmodel == 1)
    p.TMname = 'MITgcm_2.8deg';
    p.pathMatrix   = strcat(path,'/../TMs/MITgcm_2.8deg/Matrix5/TMs/matrix_nocorrection_');
    p.pathBoxes     = strcat(path,'/../TMs/MITgcm_2.8deg/Matrix5/Data/boxes.mat');
    p.pathGrid      = strcat(path,'/../TMs/MITgcm_2.8deg/grid.mat');
    p.pathConfigData = strcat(path,'/../TMs/MITgcm_2.8deg/config_data.mat');
    p.pathTemp      = strcat(path,'/../TMs/MITgcm_2.8deg/BiogeochemData/Theta_bc.mat'); 
    p.pathN0        = strcat(path,'/../TMs/MITgcm_2.8deg/N0');
    p.pathSi0       = strcat(path,'/../TMs/MITgcm_2.8deg/Si0');
    p.pathInit      = strcat(sprintf('TMs/globalInitMITgcm_%02i',length(p.u0)));
    p.pathPARday    = strcat(path,'/../TMs/MITgcm_2.8deg/parday.mat');
    
    p.dt = 0.1; % For Euler time stepping
    p.dtTransport = 0.5; % The TM time step (in units of days)
elseif nTMmodel == 2
    p.TMname = 'MITgcm_ECCO';
    p.pathMatrix = strcat(path,'/../TMs/MITgcm_ECCO/Matrix1/TMs/matrix_nocorrection_');
    p.pathBoxes = strcat(path,'/../TMs/MITgcm_ECCO/Matrix1/Data/boxes.mat');
    p.pathGrid = strcat(path,'/../TMs/MITgcm_ECCO/grid.mat');
    p.pathConfigData = strcat(path,'/../TMs/MITgcm_ECCO/config_data.mat');
    p.pathTemp = strcat(path,'/../TMs/MITgcm_ECCO/BiogeochemData/Theta_bc.mat'); 
    p.pathN0    = strcat(path,'/../TMs/MITgcm_ECCO/N0');
    p.pathSi0    = strcat(path,'/../TMs/MITgcm_ECCO/Si0');
    p.pathInit = strcat(sprintf('Transport matrix/globalInitMITgcm_ECCO_%02i',length(p.u0)));
    p.pathPARday= strcat(path,'/../TMs/MITgcm_ECCO/parday.mat'); % No parday file for high-res runs
    p.dt = 0.1; % For Euler time stepping
    p.dtTransport = 0.5; % The TM time step (in units of days)
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
    
    p.dt = 0.1; % For Euler time stepping (in units of days)
end
%
% Test that the TMs are available:
%
sTest = strcat(path,'/../TMs/',p.TMname);
if ~exist( sTest )
    error('Transport matrix directory %s does not exist.\nDownload TMs from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/',...
    sTest);
end
check(p.pathBoxes);
check(p.pathGrid);
check(p.pathConfigData);
check(p.pathTemp);
%
% Numerical parameters:
%
p.tEnd = 365; % In days
p.tSave = 365/12; % How often to save results (monthly)
p.bTransport = true; % Whether to do the transport with the transport matrix
%
% Bottom BC for nutrients:
%
p.BCmixing = [1, 0, 1]/365; % Rate of mixing nutrients into the bottom cell (1/day)
p.BCvalue = 0*p.u0 - 1; % Use the initial value concentration of the bottom concentration
p.BC_POMclosed = false; % Whether the bottom BC for POM is open or closed
%
% Light environment:
%
p.bUse_parday_light = true; % Using the parday file includes changes in cloud cover
                             % but only works with MITgcm_2.8
p.kw = 0.07; % Damping of light by water; m^-1
% Parameters used to calculate light if not using parday:
p.EinConv = 4.57; % conversion factor from W m^-2 to \mu mol s^-1 m^-2 (Thimijan & Heins 1983)
p.PARfrac = 0.4; % Fraction of light available as PAR. Source unknown

end