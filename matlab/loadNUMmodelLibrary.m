%
% Load fortran library. If it is already loaded, it returns the library
% name.
%
% In:
%  bParallel: whether to prepare the library for parallel processing
%             (default to false).
%
% Out:
%  The library name and the path of the NUMmodel library
%
function [sLibname, path] = loadNUMmodelLibrary(bParallel)

if nargin==0
    bParallel = false;
end
%
% Check matlab version:
%
if verLessThan('Matlab','9.1')
    error('Needs Matlab version 9.1 or higher.');
end
%
% Check that input files are available:
%
if ~exist('../input/input.yaml','file')
    error('The input file ../input/input.yaml is not available.');
end
%
% Find the correct library for the OS:
%
switch computer('arch')
    case {'maci','maci64','maca64'} % Note: does not distinguish between intel and arm libraries
        sLibname = 'libNUMmodel_matlab';
        sExtension = '.dylib';
    case {'glnx86','glnxa64'}
        sLibname = 'libNUMmodel_matlab';
        sExtension = '.so';
    case {'win32','win64'}
        sLibname = 'libNUMmodel_matlab';
        sExtension = '.dll';
    otherwise
        error('Architecture %s not found.\n', computer('arch'));
end
%
% Load library:
%
path = fileparts(mfilename('fullpath'));

if bParallel
    if ~exist("gcp")
        error('Running parallel simulation requires that the matlab parallel toolbox is installed.\n');
    end
    if isempty(gcp('nocreate'))
        parpool('AttachedFiles',...
            {strcat(path,'/../lib/',sLibname,sExtension),...
            strcat(path,'/../Fortran/NUMmodel_wrap_colmajor4matlab.h')});
    end
else
    if ~libisloaded(sLibname)
        [~,~] = loadlibrary(...
            strcat(path,'/../lib/',sLibname,sExtension), ...
            strcat(path,'/../Fortran/NUMmodel_wrap_colmajor4matlab.h'));
    end
end

path = regexprep(path,'(\w+)$',''); % Removes "matlab"