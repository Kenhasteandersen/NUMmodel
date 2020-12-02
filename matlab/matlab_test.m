% This uses Matlab's loadlibrary and calllib to directly call the wrapped function.  It
% does some magic with returning all the pointer variables so if they get modified they
% can be used again.
%
% This apporach is essentially the same as Julia's and python-cffi, except Matlab does the
% parsing of the h-file for you.  Sadly Octave does not yet support this at all.  Bindings
% could also be done with SWIG or through MEX-files:
% https://publicwiki.deltares.nl/display/OET/Matlab+interfacing+fortran

% Matlab reference:
% http://www.mathworks.co.uk/help/matlab/ref/loadlibrary.html#btjfvd3

libname = 'libNUMmodel_colmajor';
if libisloaded(libname)
    unloadlibrary(libname);
end
[notfound,warnings]=loadlibrary(strcat('../Fortran/',libname,'.so'), ...
    '../Fortran/NUMmodel_wrap_colmajor4matlab.h');

n = int32(10);
u = linspace(0,5,n);
dudt = 0*u;
[u,dudt] = calllib(libname, 'f_calcrates',n,u,dudt);
dudt


p = parameters;
