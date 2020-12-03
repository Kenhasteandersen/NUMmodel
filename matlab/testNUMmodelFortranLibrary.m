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

%loadLibrary();


p = parameters;
%%
setParameters(p);


u = linspace(0,10,p.n);
dudt = 0*u;
dudt = calllib(loadNUMmodelLibrary(), 'f_calcrates',10., 10.,p.n,u, 1.0, 1.0, dudt);
dudt
