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
    p.pathN0        = strcat(path,'/../TMs/MITgcm_N0');
    p.pathInit      = strcat(sprintf('Transport matrix/globalInitMITgcm_%02i',length(p.u0)));
    
    p.dt = 0.1; % For Euler time stepping
    p.dtTransport = 0.5; % The TM time step (in units of days)
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
p.bTransport = true;
%
% Set minimum concentrations:
%
p.umin = 0*p.u0;
for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    if p.typeGroups(iGroup) < 10
        p.umin(ix) = 1e-5*p.mDelta(ix(1))/p.m(ix(1)); % Minimum B concentration for unicellular groups
    end
    if p.typeGroups(iGroup)>=10 && p.typeGroups(iGroup)<100
        p.umin(ix(1)) = 1e-5*p.mDelta(ix(1))/p.m(ix(1)); % Send in some nauplii in copepod groups
    end
end
%
% Light environment:
%
p.bUse_parday_light = false; % Using the parday file includes changes in cloud cover
                             % but only works with MITgcm_2.8
% Parameters used to calculate light if not using parday:
p.EinConv = 4.57; % conversion factor from W m^-2 to \mu mol s^-1 m^-2 (Thimijan & Heins 1983)
p.PARfrac = 0.4; % Fraction of light available as PAR. Source unknown
p.kw = 0.05; % Damping of light by water; m^-1



end