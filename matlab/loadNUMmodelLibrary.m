%
% Load fortran library. If it is already loaded, it returns the library
% name.
%
function sLibname = loadNUMmodelLibrary(bParallel)

arguments
    bParallel logical = false;
end

sLibname = 'NUMmodel_matlab';

path = fileparts(mfilename('fullpath'));

if bParallel
    if isempty(gcp('nocreate'))
        parpool('AttachedFiles',...
            {strcat(path,'/../Fortran/NUMmodel_matlab.so'),...
            strcat(path,'/../Fortran/NUMmodel_wrap_colmajor4matlab.h')});
    end
else
    if ~libisloaded(sLibname)
        [notfound,warnings] = loadlibrary(...
            strcat(path,'/../Fortran/',sLibname,'.so'), ...
            strcat(path,'/../Fortran/NUMmodel_wrap_colmajor4matlab.h'));
    end
end