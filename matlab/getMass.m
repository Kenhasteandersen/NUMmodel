function m = getMass(p)

m = zeros(1,p.n);
m = calllib(loadNUMmodelLibrary(), 'f_getmass', m);
m(1:p.idxB-1) = 0;

