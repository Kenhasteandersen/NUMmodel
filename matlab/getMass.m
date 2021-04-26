function p = getMass(p)

p.m = zeros(1,p.n);
p.mDelta = p.m;
[p.m, p.mDelta] = calllib(loadNUMmodelLibrary(), 'f_getmass', p.m, p.mDelta);
p.m = logspace(-8.5, 0.1, p.n);
p.m(1:p.idxB-1) = 0;

x = log(p.m);
p.mDelta = p.m+exp(gradient(x));
