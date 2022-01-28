%
% Load fortran library. If it is already loaded, it returns the library
% name.
%
function sLibname = loadNUMmodelLibrary(bParallel)

arguments
    bParallel logical = false;
end

sLibname = 'libNUMmodel_matlab';

switch computer('arch')
    case {'maci','maci64'}
        sExtension = '.so';
    case {'glnx86','glnxa64'}
        sExtension = '.so';
    case {'win32','win64'}
        sExtension = '.dll';
end

path = fileparts(mfilename('fullpath'));

if bParallel
    if isempty(gcp('nocreate'))
        parpool('AttachedFiles',...
            {strcat(path,'/../lib/',sLibname,sExtension),...
            strcat(path,'/../Fortran/NUMmodel_wrap_colmajor4matlab.h')});
    end
else
    if ~libisloaded(sLibname)
        [notfound,warnings] = loadlibrary(...
            strcat(path,'/../lib/',sLibname,sExtension), ...
            strcat(path,'/../Fortran/NUMmodel_wrap_colmajor4matlab.h'));
    end
end