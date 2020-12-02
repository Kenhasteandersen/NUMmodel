function sLibname = loadNUMmodelLibrary

sLibname = 'libNUMmodel_colmajor';
% Should be removed when library is mature
if libisloaded(sLibname)
    unloadlibrary(sLibname);
end

[notfound,warnings] = loadlibrary(strcat('../Fortran/',sLibname,'.so'), ...
    '../Fortran/NUMmodel_wrap_colmajor4matlab.h');
