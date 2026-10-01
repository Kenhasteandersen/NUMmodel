%
% Sets the parameters for the global model simlations
%
% In:
%  p - parameter structure as returned by a setup-function (e.g. setupNUMmodel)
%  nTMmodel - Which transport matrices to use:
%       1 = MITgcm_2.8 (default); 2.8 degrees, 15 layers, 52749 wet boxes
%       2 = MITgcm_ECCO; 1 degree, 23 layers, 682604 wet boxes
%       3 = UVicOSUpicdefault; 19 layers, 87307 wet boxes. This one runs with a
%           transport time step of one day instead of half a day; see below.
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
    p.bUse_parday_light = true;

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
    p.pathPARday= strcat(path,'/../TMs/MITgcm_ECCO/parday.mat');
    p.bUse_parday_light = true;
    p.dt = 0.1; % For Euler time stepping
    p.dtTransport = 0.5; % The TM time step (in units of days)
elseif nTMmodel == 3
    p.TMname = 'UVicOSUpicdefault';
    p.pathMatrix = strcat(path,'/../TMs/UVicOSUpicdefault/Matrix1/TMs/matrix_nocorrection_');
    p.pathBoxes = strcat(path,'/../TMs/UVicOSUpicdefault/Matrix1/Data/boxes.mat');
    p.pathGrid = strcat(path,'/../TMs/UVicOSUpicdefault/grid.mat');
    p.pathConfigData = strcat(path,'/../TMs/UVicOSUpicdefault/config_data.mat');
    p.pathTemp = strcat(path,'/../TMs/UVicOSUpicdefault/BiogeochemData/Theta_bc.mat');
    % UVic comes without initial nutrient fields; simulateGlobal falls back to
    % the uniform concentrations in p.u0 when these files are not there:
    p.pathN0    = strcat(path,'/../TMs/UVicOSUpicdefault/N0');
    p.pathSi0   = strcat(path,'/../TMs/UVicOSUpicdefault/Si0');
    p.pathInit = strcat(sprintf('Transport matrix/globalInitUVicOSUpicdefault_%02i',length(p.u0)));
    % UVic does not come with a parday file; it is made by interpolating the
    % ECCO one with interpolateTMfields in "Helper functions". Use it when it is
    % there and fall back on calculating the light from the latitude and the day
    % of the year when it is not:
    p.pathPARday = strcat(path,'/../TMs/UVicOSUpicdefault/parday.mat');
    p.bUse_parday_light = isfile(p.pathPARday);

    % UVic has deltaT = 8 hours, which is the natural transport step: the
    % implicit matrix is then used as it is supplied, raised to the power one. A
    % step of half a day would need the power 1.5, which a sparse matrix cannot
    % be raised to, and a whole day needs the power 3 and moves the biomass
    % field by 11% compared to the eight-hour step.
    p.dtTransport = 1/3;
    % dtTransport has to be a whole number of Euler steps, so 1/12 of a day
    % rather than the 0.1 used with the MITgcm matrices. 1/12 gives four steps,
    % 0.1 would give three and run the biology 10% slow:
    p.dt = 1/12; % For Euler time stepping (in units of days)
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
% p.bUse_parday_light is set together with the transport matrix above, because it
% needs a parday file for the grid. Using it includes changes in cloud cover.
p.kw = 0.1; % Damping of light by water; m^-1
% Parameters used to calculate light if not using parday:
p.EinConv = 4.57; % conversion factor from W m^-2 to \mu mol s^-1 m^-2 (Thimijan & Heins 1983)
p.PARfrac = 0.4; % Fraction of light available as PAR. Source unknown

end