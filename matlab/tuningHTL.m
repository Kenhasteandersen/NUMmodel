function [y,fval,exitflag,output] = tuningHTL(parinit, siminit)
arguments
    parinit = [-1, 0.03]; % log10(mHTL) and mortHTL
    siminit = [];
end

%p = setupNUMmodel(bParallel= true);
p = setupNUMmodel([0.5 2], [10 1000], 6,6,1, bParallel= true);
p = parametersGlobal(p,1);
p.BC_POMclosed = false;
p.BCdiffusion = [1 0 10];
setSinkingPOM(p, 15)
p.tEnd = 3*365;
p.bUse_parday_light = true;
p.kw = 0.083;

mHTL = parinit(1);
mortHTL = 10^parinit(2);
bHTLdecline = false;
bHTLquadratic = true;

%
% Run a long initial condition:
%
if isempty( siminit )
    setHTL(mortHTL, mHTL, bHTLquadratic, bHTLdecline)
    p.tEnd = 10*365;
    %load('NUMmodelsmall.mat','sim');
    siminit = simulateGlobal(p,bCalcAnnualAverages=true);
    save("NUMmodel.mat",'siminit','-v7.3')
    p.tEnd = 3*365;
end

p.tSave = 1;
p.tEnd = 12;

objExpected = [1 1 1 200 800 1000];

options = optimoptions('fmincon');
options.Display = 'iter';

parmin = [-2, 0.001];
parmax = [1, 0.1];
[y,fval,exitflag,output] = fmincon(@objective, parinit, [],[],[],[], parmin,parmax, [], options);

    function err = objective(par)
        % Set parameters:
        %p.kw = par(1);
        %setSinkingPOM(p,par(2));
        setHTL(10^par(1), par(2), bHTLquadratic, bHTLdecline);

        % Simulate:
        sim = simulateGlobal(p, siminit,bCalcAnnualAverages=true);%, bVerbose=false);

        % Plots:
        figure(1)
        obj = EvaluateRun(sim);
        figure(2)%f
        plotWatercolumnTime(sim,60,-15,depthMax=200)
        drawnow

        % Evaluate objective criterion:
        obj = obj(1:3) ./ objExpected(1:3); % Evaluate only using biomasses and not NPP
        err = double(mean(abs(log(obj))));
        fprintf("Pars: %f, %f; err: %f\n",[par,err]);
        fprintf("Objective vector: %f, %f, %f, %f, %f, %f\n", obj)
        fprintf("====================\n");
    end

end