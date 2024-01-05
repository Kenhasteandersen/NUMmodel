p = parametersGlobal( setupNUMmodel([0.5 2], [10 1000], 6,6,1, bParallel= true) );

%p.tEnd = 12;
%p.tSave = 1;
maxDepth = [100,200,300,500,1000,5000];

sim = {};
for i = 1:length(maxDepth)
    p.maxDepth = maxDepth(i);

    tic
    sim{i} = simulateGlobal(p, bVerbose=false, bCalcAnnualAverages=true);
    time(i) = toc
    status(i,:) = EvaluateRun(sim{i});
end
