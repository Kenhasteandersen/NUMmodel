function sLibname = loadNUMmodelLibrary

sLibname = 'libNUMmodel_colmajor';
if ~libisloaded(sLibname)
    [notfound,warnings] = loadlibrary(strcat('../Fortran/',sLibname,'.so'), ...
        '../Fortran/NUMmodel_wrap_colmajor4matlab.h');
end
