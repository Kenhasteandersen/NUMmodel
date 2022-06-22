p = setupGeneralistsDiatoms;
p = parametersChemostat(p);
p.tEnd = 365;

    setHTL(0.1, 1/500^1.5,false,false);

d=logspace(-4,-1,10);
L=linspace(0,500,11);

tic
for i=1:length(d)
simd(i) = simulateChemostat(p, d(i));
end
toc