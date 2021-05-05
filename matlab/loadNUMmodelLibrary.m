%
% Load fortran library. If it is already loaded, it returns the library
% name.
%
function sLibname = loadNUMmodelLibrary(bParallel)

arguments
    bParallel logical = false;
end

sLibname = 'NUMmodel_matlab';

if bParallel
    if isempty(gcp('nocreate'))
        parpool('AttachedFiles',...
            {'../Fortran/NUMmodel_matlab.so',...
            '../Fortran/NUMmodel_wrap_colmajor4matlab.h'});
    end
else
    if ~libisloaded(sLibname)
        [notfound,warnings] = loadlibrary(strcat(...
            '../Fortran/',sLibname,'.so'), ...
            '../Fortran/NUMmodel_wrap_colmajor4matlab.h');
    end
end