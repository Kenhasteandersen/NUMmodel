function sLibname = loadNUMmodelLibrary

sLibname = 'NUMmodel';
if ~libisloaded(sLibname)
    [notfound,warnings] = loadlibrary(strcat('../Fortran/',sLibname,'.so'), ...
        '../Fortran/NUMmodel_wrap_colmajor4matlab.h');
    calllib(loadNUMmodelLibrary(), 'f_parametersgeneralistsonly');
end
